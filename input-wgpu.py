###############################################################################
# Some standard modules
import os
import sys
import copy
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
###############################################################################

import openmc
import json
import argparse



import helper



cspath = "/openmc/openmc-data/v0.12/lib80x_hdf5/cross_sections.xml"
basepath = "fetter-20210728"
particles = 1
batches = 100
ages = [0, 5, 10, 15, 20, 25, 30, 35, 40]
ages = [float(a) for a in ages]

parser = argparse.ArgumentParser()
parser.add_argument("-p", "--plot", action="store_true",
                    help="plot some cuts through the geometry")
parser.add_argument("-b", "--basepath", 
                    help="set basepath for calculation")
parser.add_argument("-c", "--batches", type = int,
                    help="set number of batches")
parser.add_argument("-n", "--particles", type = int,
                    help="set number of particles")
parser.add_argument("-a", "--ages",
                    help="set ages for weapon in calculation")

args = parser.parse_args()
if not args.basepath is None:
    basepath = args.basepath
if not args.ages is None:
    stringages = args.ages.split()
    ages = [float(a) for a in stringages]
if not args.batches is None:
    batches = args.batches
if not args.particles is None:
    particles = args.particles
if args.plot:
    plot = True

geometries = ["pit-outer",
              "ref-outer",
              "tam-outer",
              "exp-outer",
              "full"]

evpath = os.path.join(basepath, "eigenvalue")
fspath = os.path.join(basepath, "fixed-source")
if not os.path.exists(basepath):
    os.mkdir(basepath)
    os.mkdir(evpath)
    os.mkdir(fspath)
    for geo in geometries:
        os.mkdir(os.path.join(evpath, geo))
        os.mkdir(os.path.join(fspath, geo))
    for age in ages:
        os.mkdir(os.path.join(evpath, "full-{:.2f}".format(age)))
        os.mkdir(os.path.join(fspath, "full-{:.2f}".format(age)))

# dump basic settings to json file
with open(os.path.join(basepath, "calculation.json"), 'w') as f:
    json.dump({'ages': ages,
               'geometries': geometries,
               'particles': particles,
               'basepath': basepath,
               'cspath': cspath,
               'batches': batches},
              f)
    f.close()
    
###############################################################################
# Materials
###############################################################################

# need pit radius already here
pitOR = 5
pitThickness = 0.75
pitIR = pitOR - pitThickness

puvec = helper.puvector()
puvec.setwo('Pu238', 0.00005)
puvec.setwo('Pu239', 0.993)
puvec.setwo('Pu240', 0.060)
puvec.setwo('Pu241', 0.0044)
puvec.setwo('Pu242', 0.00015)
puvec.setwo('Am241', 0)
puvec.createagedvector(ages)

pudf = puvec.pudf
puagedf = puvec.puagedf

wpuRho = 4000 / (4/3 * np.pi * (pitOR ** 3 - pitIR ** 3)) # density based on Fetter Model volume / mass
wpuMat = {}
for age in ages:
    wpuMat[age] = openmc.Material(name = "Weapon-grade Plutonium - {:.2f} years old".format(age))
    wpuMat[age].set_density("g/cm3", wpuRho)
    wpuMat[age].add_nuclide("Pu238", puagedf.loc[('Pu238', age), 'wo'], "wo")
    wpuMat[age].add_nuclide("Pu239", puagedf.loc[('Pu239', age), 'wo'], "wo")
    wpuMat[age].add_nuclide("Pu240", puagedf.loc[('Pu240', age), 'wo'], "wo") #should be wo all along
    wpuMat[age].add_nuclide("Pu241", puagedf.loc[('Pu241', age), 'wo'], "wo") # need to consider for decay!
    wpuMat[age].add_nuclide("Pu242", puagedf.loc[('Pu242', age), 'wo'], "wo")
    wpuMat[age].add_nuclide("Am241", puagedf.loc[('Am241', age), 'wo'], "wo")
    wpuMat[age].add_element("O", 0.002, "wo")

berMat = openmc.Material(name = "Beryllium Reflector")
berMat.set_density("g/cm3", 1.848) #PNNL Compendium page 37
berMat.add_element("Be", 1)

tunMat = openmc.Material(name = "Tungsten Tamper")
tunMat.set_density("g/cm3", 19.3) #PNNL Compendium page 318
tunMat.add_element("W", 1)

hmxMat = openmc.Material(name = "HMX explosive")
hmxMat.set_density("g/cm3", 1.890) #PNNL Compendium page 125
hmxMat.add_element("H", 0.285714, "ao")
hmxMat.add_element("C", 0.142857, "ao")
hmxMat.add_element("N", 0.285714, "ao")
hmxMat.add_element("O", 0.285714, "ao")

aluMat = openmc.Material(name = "Aluminum Case")
aluMat.set_density("g/cm3", 2.6989) #PNNL Compendium page 20
aluMat.add_element("Al", 1)

airMat = openmc.Material(name = "Air")
airMat.set_density("g/cm3", 0.001205) #PNNL Compendium page 18
airMat.add_element("C", 0.000150, "ao")
airMat.add_element("N", 0.784431, "ao")
airMat.add_element("O", 0.210748, "ao")
airMat.add_element("Ar", 0.004671, "ao")

materiallist = [wpuMat[0], # default age = 0 years
                berMat, tunMat, hmxMat, aluMat, airMat]
materials = openmc.Materials(materiallist)
materials.cross_sections = cspath

for geo in geometries:
    materials.export_to_xml(os.path.join(evpath, geo, "materials.xml"))
    materials.export_to_xml(os.path.join(fspath, geo, "materials.xml"))

for age in ages:
    materiallist = [wpuMat[age], 
                    berMat, tunMat, hmxMat, aluMat, airMat]
    materials = openmc.Materials(materiallist)
    materials.cross_sections = cspath
    materials.export_to_xml(os.path.join(evpath, "full-{:.2f}".format(age), "materials.xml"))
    materials.export_to_xml(os.path.join(fspath, "full-{:.2f}".format(age), "materials.xml"))


###############################################################################
# Geometry
###############################################################################

# Pit

cenSurface = openmc.Sphere(r = pitIR)
cenCell = openmc.Cell()
cenCell.region = -cenSurface
cenCell.name = "Center of weapon"

pitSurface = openmc.Sphere(r = pitOR)
pitCell = openmc.Cell()
pitCell.region = +cenSurface & -pitSurface
pitCell.name = "Pit"

# Reflector
refThickness = 2
refOR = pitOR + refThickness

refSurface = openmc.Sphere(r = refOR)
refCell = openmc.Cell()
refCell.region = +pitSurface & -refSurface
refCell.name = "Reflector"

# Tamper
tamThickness = 3
tamOR = refOR + tamThickness

tamSurface = openmc.Sphere(r = tamOR)
tamCell = openmc.Cell()
tamCell.region = +refSurface & -tamSurface
tamCell.name = "Tamper"

# Explosive
expThickness = 10
expOR = tamOR + expThickness

expSurface = openmc.Sphere(r = expOR)
expCell = openmc.Cell()
expCell.region = +tamSurface & -expSurface
expCell.name = "Conventional Explosive"

# Aluminum Case
aluThickness = 1
aluOR = expOR + aluThickness

aluSurface = openmc.Sphere(r = aluOR)
aluCell = openmc.Cell()
aluCell.region = +expSurface & -aluSurface
aluCell.name = "Aluminum Case"

# Outside
outThickness = 20
outOR = aluOR + outThickness

outSurface = openmc.Sphere(r = outOR, boundary_type='vacuum')
outCell = openmc.Cell()
outCell.region = +aluSurface & -outSurface
outCell.name = "Outside (vacuum)"

# Create multiple geometries
pitCell.fill = wpuMat[0] 
cells = [cenCell, pitCell, refCell, tamCell, expCell, aluCell, outCell]
root = openmc.Universe(cells = cells)
geometry = openmc.Geometry(root)
geometry.export_to_xml(os.path.join(evpath, "pit-outer", "geometry.xml"))
geometry.export_to_xml(os.path.join(fspath, "pit-outer", "geometry.xml"))

refCell.fill = berMat
cells = [cenCell, pitCell, refCell, tamCell, expCell, aluCell, outCell]
root = openmc.Universe(cells = cells)
geometry = openmc.Geometry(root)
geometry.export_to_xml(os.path.join(evpath, "ref-outer", "geometry.xml"))
geometry.export_to_xml(os.path.join(fspath, "ref-outer", "geometry.xml"))

tamCell.fill = tunMat
cells = [cenCell, pitCell, refCell, tamCell, expCell, aluCell, outCell]
root = openmc.Universe(cells = cells)
geometry = openmc.Geometry(root)
geometry.export_to_xml(os.path.join(evpath, "tam-outer", "geometry.xml"))
geometry.export_to_xml(os.path.join(fspath, "tam-outer", "geometry.xml"))

expCell.fill = hmxMat
cells = [cenCell, pitCell, refCell, tamCell, expCell, aluCell, outCell]
root = openmc.Universe(cells = cells)
geometry = openmc.Geometry(root)
geometry.export_to_xml(os.path.join(evpath, "exp-outer", "geometry.xml"))
geometry.export_to_xml(os.path.join(fspath, "exp-outer", "geometry.xml"))

aluCell.fill = aluMat
cells = [cenCell, pitCell, refCell, tamCell, expCell, aluCell, outCell]
root = openmc.Universe(cells = cells)
geometry = openmc.Geometry(root)
geometry.export_to_xml(os.path.join(evpath, "full", "geometry.xml"))
geometry.export_to_xml(os.path.join(fspath, "full", "geometry.xml"))

for age in ages:
    pitCell.fill = wpuMat[age]
    cells = [cenCell, pitCell, refCell, tamCell, expCell, aluCell, outCell]
    root = openmc.Universe(cells = cells)
    geometry = openmc.Geometry(root)
    geometry.export_to_xml(os.path.join(evpath, "full-{:.2f}".format(age), "geometry.xml"))
    geometry.export_to_xml(os.path.join(fspath, "full-{:.2f}".format(age), "geometry.xml"))

plt.figure(figsize=(5, 5))
root.plot(width=(47, 47))
plt.title("Geometry of Fetter model in OpenMC")
plt.savefig(os.path.join(basepath, "fetter_geometry.png"))


###############################################################################
# Tallies
###############################################################################

tallies = openmc.Tallies()

allcellfilter = openmc.CellFilter([cenCell, pitCell, refCell, tamCell, expCell, aluCell, outCell])
allcellfromfilter = openmc.CellFromFilter([cenCell, pitCell, refCell, tamCell, expCell, aluCell, outCell])
allsurfacefilter = openmc.SurfaceFilter([cenSurface, pitSurface, refSurface, tamSurface, expSurface, aluSurface, outSurface])
energyfilter = openmc.EnergyFilter(np.logspace(-3, 8, num=100))
cellbornfilter = openmc.CellbornFilter([cenCell, pitCell, refCell, tamCell, expCell, aluCell, outCell])

# Flux
tally = openmc.Tally(name='flux')
tally.filters = [allcellfilter,
                 energyfilter]
tally.scores = ['flux']
tallies.append(tally)


# Leakage with Mesh
# Instantiate a tally mesh
mesh = openmc.RegularMesh(mesh_id=1)
mesh.dimension = [1, 1, 1]
mhw = aluOR + 0.01            
mesh.lower_left = [-mhw, -mhw, -mhw]
mesh.width = [2 * mhw, 2 * mhw, 2 * mhw]
meshsurface_filter = openmc.MeshSurfaceFilter(mesh)

# Instantiate thermal, fast, and total leakage tallies
leak = openmc.Tally(name='leakage')
leak.filters = [meshsurface_filter]
leak.scores = ['current']
tallies.append(leak)

thermal_leak = openmc.Tally(name='thermal leakage')
thermal_leak.filters = [meshsurface_filter, openmc.EnergyFilter([0., 0.625])]
thermal_leak.scores = ['current']
tallies.append(thermal_leak)

fast_leak = openmc.Tally(name='fast leakage')
fast_leak.filters = [meshsurface_filter, openmc.EnergyFilter([0.625, 20.0e6])]
fast_leak.scores = ['current']
tallies.append(fast_leak)

# K-Eigenvalue (infinity) tallies
fiss_rate = openmc.Tally(name='fiss. rate')
abs_rate = openmc.Tally(name='abs. rate')
fiss_rate.scores = ['nu-fission']
abs_rate.scores = ['absorption']
tallies += (fiss_rate, abs_rate)

nuscatter_rate = openmc.Tally(name='nu scatter')
nuscatter_rate.scores = ['nu-scatter']
tallies.append(nuscatter_rate)

rates = openmc.Tally(name='various reaction rates')
rates.scores = ['absorption', 'fission', 'scatter', 'total', 
                '(n,2nd)', '(n,2n)', '(n,3n)', '(n,na)', '(n,n3a)', '(n,2na)',
                '(n,3na)', '(n,np)', '(n,n2a)', '(n,2n2a)', '(n,nd)', '(n,nt)',
                '(n,n3He)', '(n,nd2a)', '(n,nt2a)', '(n,4n)', '(n,2np)',
                '(n,3np)', '(n,n2p)',
                '(n,gamma)', '(n,p)', '(n,t)', '(n,3He)', '(n,a)',
                '(n,2a)', '(n,3a)', '(n,2p)', '(n,pa)', '(n,t2a)', '(n,d2a)',
                '(n,pd)', '(n,pt)', '(n,da)',
                'nu-fission', 'nu-scatter'
                ]
tallies.append(rates)

# Resonance Escape Probability tallies
therm_abs_rate = openmc.Tally(name='therm. abs. rate')
therm_abs_rate.scores = ['absorption']
therm_abs_rate.filters = [openmc.EnergyFilter([0., 0.625])]
tallies.append(therm_abs_rate)

# Thermal Flux Utilization tallies
fuel_therm_abs_rate = openmc.Tally(name='fuel therm. abs. rate')
fuel_therm_abs_rate.scores = ['absorption']
fuel_therm_abs_rate.filters = [openmc.EnergyFilter([0., 0.625]),
                               openmc.CellFilter([pitCell])]
tallies.append(fuel_therm_abs_rate)

# Fast Fission Factor tallies
therm_fiss_rate = openmc.Tally(name='therm. fiss. rate')
therm_fiss_rate.scores = ['nu-fission']
therm_fiss_rate.filters = [openmc.EnergyFilter([0., 0.625])]
tallies.append(therm_fiss_rate)

# Current through surfaces
tally = openmc.Tally(name='surface current')
tally.filters = [allsurfacefilter, energyfilter, allcellfromfilter]
tally.scores = ['current']
tallies.append(tally)

# Current with cells, fromcells and surface filters
tally = openmc.Tally(name='plain surface current')
tally.filters = [allsurfacefilter]
tally.scores = ['current']
tallies.append(tally)

tally = openmc.Tally(name='plain cell current')
tally.filters = [allcellfilter]
tally.scores = ['current']
tallies.append(tally)

tally = openmc.Tally(name='plain from cell current')
tally.filters = [allcellfromfilter]
tally.scores = ['current']
tallies.append(tally)

tally = openmc.Tally(name='leaked neutrons with birth location')
tally.filters = [allcellfilter, cellbornfilter]
tally.scores = ['flux']
tallies.append(tally)

tally = openmc.Tally(name='total flux')
tally.scores = ['flux']
tallies.append(tally)

tally = openmc.Tally(name='cellwise flux')
tally.filters = [allcellfilter]
tally.scores = ['flux']
tallies.append(tally)

for geo in geometries:
    tallies.export_to_xml(os.path.join(fspath, geo, "tallies.xml"))
    tallies.export_to_xml(os.path.join(evpath, geo, "tallies.xml"))
for age in ages:
    tallies.export_to_xml(os.path.join(fspath, "full-{:.2f}".format(age)))
    tallies.export_to_xml(os.path.join(evpath, "full-{:.2f}".format(age)))
                 
###############################################################################
# Source & Settings
###############################################################################


bounds = [-pitOR, -pitOR, -pitOR, pitOR, pitOR, pitOR]
uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable=True)
evsrc = openmc.Source(space = uniform_dist, particle = 'neutron')


settings = openmc.Settings()
settings.run_mode = 'eigenvalue'
settings.batches = batches + 20
settings.inactive = 20
settings.particles = particles

settings.source = evsrc
for geo in geometries:
    settings.export_to_xml(os.path.join(evpath, geo, "settings.xml"))
for age in ages:
    settings.export_to_xml(os.path.join(evpath, "full-{:.2f}".format(age)))

settings.run_mode = 'fixed source'
settings.inactive = 0
settings.batches = batches
fssrc = []

normalizationconstant = sum(pudf['sfneutrons'] * pudf['wo'])
for iso in pudf.index:
    tmpsrc = openmc.Source(space = uniform_dist, particle = 'neutron')
    tmpsrc.energy = openmc.stats.Watt(a=1/pudf.loc[iso, 'watt-a'], b=pudf.loc[iso, 'watt-b'])
    tmpsrc.strength = pudf.loc[iso, 'sfneutrons'] * pudf.loc[iso, 'wo'] / normalizationconstant
    fssrc.append(tmpsrc)
    
settings.source = fssrc
for geo in geometries:
    settings.export_to_xml(os.path.join(fspath, geo, "settings.xml"))

for age in ages:
    fssrc = []
    subset = puagedf.loc[(slice(None), age), :]
    normalizationconstant = sum(subset['wo'] * subset['sfneutrons'])
    for iso in pudf.index:
        tmpsrc = openmc.Source(space = uniform_dist, particle = 'neutron')
        tmpsrc.energy = openmc.stats.Watt(a=1/pudf.loc[iso, 'watt-a'], b=pudf.loc[iso, 'watt-b'])
        tmpsrc.strength = puagedf.loc[(iso, age), 'sfneutrons'] * puagedf.loc[(iso, age), 'wo'] / normalizationconstant
        fssrc.append(tmpsrc)
        settings.source = fssrc
    
        settings.export_to_xml(os.path.join(fspath, "full-{:.2f}".format(age)))

mass = 4000        
for age in ages:
    sourcesum = 0
    for iso in pudf.index:
        sourcesum += puagedf.loc[(iso, age), 'sfneutrons'] * puagedf.loc[(iso, age), 'wo']
    print("A 4kg WPu weapon ({:.2f} years old) emits {:f} neutrons/s".format(age, sourcesum * 4000))

