import json
import os
import matplotlib.pyplot as plt
import math
import numpy as np
import sys
import argparse
import shutil

import openmc
import helper

cspath = "/openmc/openmc-data/v0.12/lib80x_hdf5/cross_sections.xml"
cspathfetter = cspath
basepath = "simple-cosmic-ray-neutrons"
particles = 10000
batches = 10
plot = False
survival = False
discard = False
maxenergy = 2e7
weaponage = 0
weightcutoff = 0.25

parser = argparse.ArgumentParser()
parser.add_argument("-p", "--plot", action="store_true",
                    help="plot some cuts through the geometry")
parser.add_argument("-s", "--survival", action="store_true",
                    help="turn on survival biasing")
parser.add_argument("-b", "--basepath", 
                    help="set basepath for calculation")
parser.add_argument("-c", "--batches", type = int,
                    help="set number of batches")
parser.add_argument("-m", "--maxenergy", type = float,
                    help="source particles cutoff energy?")
parser.add_argument("-d", "--discard", action="store_true",
                    help="discard source particles above cutoff energy?")
parser.add_argument("-x", "--xsection",
                    help="specify cross section path for cosmic ray simulations")
parser.add_argument("-y", "--xsectionfetter",
                    help="specify cross section path for fetter source simulations")
parser.add_argument("-w", "--weightcutoff", type = float,
                    help="weight cutoff value (only used if survial biasing is on)")

args = parser.parse_args()
if args.plot:
    plot = True
if args.survival:
    survival = True
if not args.weightcutoff is None:
    weightcutoff = args.weightcutoff
if not args.basepath is None:
    basepath = args.basepath
if not args.batches is None:
    batches = args.batches
if not args.maxenergy is None:
    maxenergy = args.maxenergy
if args.discard:
    discard = True
if not args.xsection is None:
    cspath = args.xsection
if not args.xsectionfetter is None:
    cspathfetter = args.xsectionfetter

concretepath = os.path.join(basepath, "concrete")
waterpath = os.path.join(basepath, "water")
soilpath = os.path.join(basepath, "soil")
concrete2100path = os.path.join(basepath, "concrete-2100")
water2100path = os.path.join(basepath, "water-2100")
soil2100path = os.path.join(basepath, "soil-2100")
concretefetterpath = os.path.join(basepath, "concrete-fetter")
waterfetterpath = os.path.join(basepath, "water-fetter")
soilfetterpath = os.path.join(basepath, "soil-fetter")
if not os.path.exists(basepath):
    os.mkdir(basepath)
    os.mkdir(concretepath)
    os.mkdir(waterpath)
    os.mkdir(soilpath)
    os.mkdir(concrete2100path)
    os.mkdir(water2100path)
    os.mkdir(soil2100path)
    os.mkdir(concretefetterpath)
    os.mkdir(waterfetterpath)
    os.mkdir(soilfetterpath)

# dump basic settings to json file
with open(os.path.join(basepath, "calculation.json"), 'w') as f:
    json.dump({'particles': particles,
               'basepath': basepath,
               'cspath': cspath,
               'cspathfetter': cspathfetter,
               'batches': batches,
               'discard': discard,
               'maxenergy': maxenergy,
               'weaponage': weaponage},
              f)
    f.close()

###############################################################################
# Materials
###############################################################################

airMat = openmc.Material(name='Air, Dry (near sea level)')
airMat.set_density('g/cm3', 1.205e-3) # PNNL Compendium page 18
airMat.add_element('C', 0.000124, 'wo')
airMat.add_element('N', 0.755268, 'wo')
airMat.add_element('O', 0.231781, 'wo')
airMat.add_element('Ar', 0.012827, 'wo')

soilMat = openmc.Material(name='Earth, Typical Western US')
soilMat.set_density('g/cm3', 1.52) # PNNL Compendium p.118 (references Brewer 2009, which has a density of 5.5g/cm3)
soilMat.add_element('H', 0.316855, 'ao')
soilMat.add_element('O', 0.501581, 'ao')
soilMat.add_element('Al', 0.039951, 'ao')
soilMat.add_element('Si', 0.141613, 'ao')

concreteRegularMat = openmc.Material(name = "Regular Concrete")
concreteRegularMat.set_density("g/cm3", 2.3) # PNNL Compendium page 112
concreteRegularMat.add_element("H", 0.16038, "ao")
concreteRegularMat.add_element("O", 0.563183, "ao")
concreteRegularMat.add_element("Na", 0.021365, "ao")
concreteRegularMat.add_element("Al", 0.021343, "ao")
concreteRegularMat.add_element("Si", 0.203231, "ao")
concreteRegularMat.add_element("Ca", 0.018595, "ao")
concreteRegularMat.add_element("Fe", 0.004246, "ao")

waterMat = openmc.Material(name = "Water")
waterMat.set_density("g/cm3", 0.998207) #PNNL Compendium page 338
waterMat.add_element("H", 2, "ao")
waterMat.add_element("O", 1, "ao")
waterMat.add_s_alpha_beta('c_H_in_H2O')

#*******************************************************************************
# fetter model

# need pit radius already here
pitOR = 5
pitThickness = 0.75
pitIR = pitOR - pitThickness

puvec = helper.createwgpuvec() # helper has the fetter model wgpu settings
ages = [weaponage]
puvec.createagedvector(ages)
puagedf = puvec.puagedf
pudf = puvec.pudf

wpuRho = 4000 / (4/3 * np.pi * (pitOR ** 3 - pitIR ** 3)) # density based on Fetter Model volume / mass
wpuMat = {}
for age in ages:
    wpuMat[age] = openmc.Material(name = "Weapon-grade Plutonium - {:02d} years old".format(age))
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

materiallist = [airMat, concreteRegularMat]
materials = openmc.Materials(materiallist)
helper.checkcrosssections(cspath, materials)
materials.cross_sections = cspath
materials.export_to_xml(os.path.join(concretepath, "materials.xml"))
materials.export_to_xml(os.path.join(concrete2100path, "materials.xml"))

materiallist = [airMat, concreteRegularMat, wpuMat[weaponage], berMat, tunMat, hmxMat, aluMat]
materials = openmc.Materials(materiallist)
helper.checkcrosssections(cspathfetter, materials)
materials.cross_sections = cspathfetter
materials.export_to_xml(os.path.join(concretefetterpath, "materials.xml"))

materiallist = [airMat, waterMat]
materials = openmc.Materials(materiallist)
helper.checkcrosssections(cspath, materials)
materials.cross_sections = cspath
materials.export_to_xml(os.path.join(waterpath, "materials.xml"))
materials.export_to_xml(os.path.join(water2100path, "materials.xml"))

materiallist = [airMat, waterMat, wpuMat[weaponage], berMat, tunMat, hmxMat, aluMat]
materials = openmc.Materials(materiallist)
helper.checkcrosssections(cspathfetter, materials)
materials.cross_sections = cspathfetter
materials.export_to_xml(os.path.join(waterfetterpath, "materials.xml"))

materiallist = [airMat, soilMat]
materials = openmc.Materials(materiallist)
helper.checkcrosssections(cspath, materials)
materials.cross_sections = cspath
materials.export_to_xml(os.path.join(soilpath, "materials.xml"))
materials.export_to_xml(os.path.join(soil2100path, "materials.xml"))

materiallist = [airMat, soilMat, wpuMat[weaponage], berMat, tunMat, hmxMat, aluMat]
materials = openmc.Materials(materiallist)
helper.checkcrosssections(cspathfetter, materials)
materials.cross_sections = cspathfetter
materials.export_to_xml(os.path.join(soilfetterpath, "materials.xml"))


###############################################################################
# Geometry
###############################################################################

hw = 10000
zmin = -500
zmax = 1500

groundBottom = openmc.ZPlane(zmin)
groundAirSurface = openmc.ZPlane(0)
airTop = openmc.ZPlane(zmax)

left = openmc.XPlane(-hw)
right = openmc.XPlane(hw)
front = openmc.YPlane(-hw)
back = openmc.YPlane(hw)

concreteCell = openmc.Cell()
concreteCell.region = +groundBottom & -groundAirSurface & +left & -right & +front & -back
concreteCell.fill = concreteRegularMat

waterCell = openmc.Cell()
waterCell.region = +groundBottom & -groundAirSurface & +left & -right & +front & -back
waterCell.fill = waterMat

soilCell = openmc.Cell()
soilCell.region = concreteCell.region = +groundBottom & -groundAirSurface & +left & -right & +front & -back
soilCell.fill = soilMat

airCell = openmc.Cell()
airCell.region = +groundAirSurface & -airTop & +left & -right & +front & -back
airCell.fill = airMat

groundBottom.boundary_type = 'vacuum'
airTop.boundary_type = 'vacuum'

left.boundary_type = 'periodic'
left.periodic_surface = right
right.boundary_type = 'periodic'
right.periodic_surface = left
front.boundary_type = 'periodic'
front.periodic_surface = back
back.boundary_type = 'periodic'
back.periodic_surface = front

root = openmc.Universe(cells = [airCell, soilCell])
geometry = openmc.Geometry(root)
geometry.export_to_xml(os.path.join(soilpath, "geometry.xml"))
geometry.export_to_xml(os.path.join(soil2100path, "geometry.xml"))

root = openmc.Universe(cells = [airCell, concreteCell])
geometry = openmc.Geometry(root)
geometry.export_to_xml(os.path.join(concretepath, "geometry.xml"))
geometry.export_to_xml(os.path.join(concrete2100path, "geometry.xml"))

root = openmc.Universe(cells = [airCell, waterCell])
geometry = openmc.Geometry(root)
geometry.export_to_xml(os.path.join(waterpath, "geometry.xml"))
geometry.export_to_xml(os.path.join(water2100path, "geometry.xml"))




#*******************************************************************************
# fetter model
# always at origin, z=1m

left.boundary_type = 'vacuum'
right.boundary_type = 'vacuum'
front.boundary_type = 'vacuum'
back.boundary_type = 'vacuum'

weaponz = 150
# Pit
cenSurface = openmc.Sphere(r = pitIR, z0=weaponz)
cenCell = openmc.Cell()
cenCell.region = -cenSurface
cenCell.name = "Center of weapon"

pitSurface = openmc.Sphere(r = pitOR, z0=weaponz)
pitCell = openmc.Cell()
pitCell.region = +cenSurface & -pitSurface
pitCell.name = "Pit"

# Reflector
refThickness = 2
refOR = pitOR + refThickness

refSurface = openmc.Sphere(r = refOR, z0=weaponz)
refCell = openmc.Cell()
refCell.region = +pitSurface & -refSurface
refCell.name = "Reflector"

# Tamper
tamThickness = 3
tamOR = refOR + tamThickness

tamSurface = openmc.Sphere(r = tamOR, z0=weaponz)
tamCell = openmc.Cell()
tamCell.region = +refSurface & -tamSurface
tamCell.name = "Tamper"

# Explosive
expThickness = 10
expOR = tamOR + expThickness

expSurface = openmc.Sphere(r = expOR, z0=weaponz)
expCell = openmc.Cell()
expCell.region = +tamSurface & -expSurface
expCell.name = "Conventional Explosive"

# Aluminum Case
aluThickness = 1
aluOR = expOR + aluThickness

aluSurface = openmc.Sphere(r = aluOR, z0=weaponz)
aluCell = openmc.Cell()
aluCell.region = +expSurface & -aluSurface
aluCell.name = "Aluminum Case"

pitCell.fill = wpuMat[weaponage] 
refCell.fill = berMat
tamCell.fill = tunMat
expCell.fill = hmxMat
aluCell.fill = aluMat

#reduce air cell by weapon
airCell.region &= +aluSurface

root = openmc.Universe(cells = [airCell, concreteCell, cenCell, pitCell, refCell, tamCell, expCell, aluCell])
geometry = openmc.Geometry(root)
geometry.export_to_xml(os.path.join(concretefetterpath, "geometry.xml"))

root = openmc.Universe(cells = [airCell, waterCell, cenCell, pitCell, refCell, tamCell, expCell, aluCell])
geometry = openmc.Geometry(root)
geometry.export_to_xml(os.path.join(waterfetterpath, "geometry.xml"))

root = openmc.Universe(cells = [airCell, soilCell, cenCell, pitCell, refCell, tamCell, expCell, aluCell])
geometry = openmc.Geometry(root)
geometry.export_to_xml(os.path.join(soilfetterpath, "geometry.xml"))

###############################################################################
# Source & Settings
###############################################################################
settings = openmc.Settings()
settings.run_mode = 'fixed source'
settings.inactive = 0
settings.batches = batches
settings.particles = particles
settings.output = {'tallies': False}
if survival:
    settings.survival_biasing = True
    cutoff = {}
    cutoff['weight'] = weightcutoff
    settings.cutoff = cutoff

#*******************************************************************************
#cosmic ray
scriptdir = os.path.dirname(os.path.realpath(__file__))
for p in [concretepath, concrete2100path, waterpath, water2100path, soilpath, soil2100path]:
    if os.path.islink(os.path.join(p, "libsource.so")):
        os.unlink(os.path.join(p, "libsource.so"))
    os.symlink(os.path.join(scriptdir, "../cry-with-openmc/build/libsource.so"), os.path.join(p, "libsource.so"))
    if os.path.islink(os.path.join(p, "data")):
        os.unlink(os.path.join(p, "data"))
    os.symlink(os.path.join(scriptdir, "../cry-with-openmc/build/CRY-1.7-prefix/src/CRY-1.7/data"), os.path.join(p, "data"))
source = openmc.Source(library = "./libsource.so")
if discard:
    discardstring = "discard"
else:
    discardstring = "limit"
latitude = 45
source.parameters = "{:f} {:s} {:f} {:f} {:f} returnNeutrons 1 returnProtons 0 returnGammas 0 returnMuons 0 returnElectrons 0 returnPions 0 date 1-1-2008 latitude {:f} altitude 0 subboxLength {:f}".format(maxenergy, discardstring, 0, 0, zmax - 0.001, latitude, 2 * hw / 100)
    
settings.source = source

for p in [concretepath, waterpath, soilpath]:
    settings.export_to_xml(os.path.join(p, "settings.xml"))

settings.source.parameters = "{:f} {:s} {:f} {:f} {:f} returnNeutrons 1 returnProtons 0 returnGammas 0 returnMuons 0 returnElectrons 0 returnPions 0 date 1-1-2008 latitude {:f} altitude 2100 subboxLength {:f}".format(maxenergy, discardstring, 0, 0, zmax - 0.001, latitude, 2 * hw / 100)

for p in [concrete2100path, water2100path, soil2100path]:
    settings.export_to_xml(os.path.join(p, "settings.xml"))

#*******************************************************************************
#fetter
sourcex = 0
sourcey = 0
sourcez = weaponz

bounds = [sourcex - pitOR, sourcey - pitOR, sourcez - pitOR,
          sourcex + pitOR, sourcey + pitOR, sourcez + pitOR]
uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable=True)
subset = puagedf.loc[(slice(None), weaponage), :]
normalizationconstant = sum(subset['wo'] * subset['sfneutrons'])
source = []
for iso in pudf.index:
    tmpsrc = openmc.Source(space = uniform_dist, particle = 'neutron')
    tmpsrc.energy = openmc.stats.Watt(a=1/pudf.loc[iso, 'watt-a'], b=pudf.loc[iso, 'watt-b'])
    tmpsrc.strength = puagedf.loc[(iso, age), 'sfneutrons'] * puagedf.loc[(iso, age), 'wo'] / normalizationconstant
    source.append(tmpsrc)
settings.source = source

for p in [concretefetterpath, waterfetterpath, soilfetterpath]:
    settings.export_to_xml(os.path.join(p, "settings.xml"))

###############################################################################
# Tallies
###############################################################################

tallies = openmc.Tallies()

xmin = -hw
xmax = +hw
ymin = -hw
ymax = +hw

xstep = 100
tallyxmin = math.floor(xmin / xstep) * xstep + xstep / 2
tallyxmax = math.ceil(xmax / xstep) * xstep - xstep / 2
ystep = 100
tallyymin = math.floor(ymin / ystep) * ystep + xstep / 2
tallyymax = math.ceil(ymax / ystep) * ystep - xstep / 2
zstep = 100
tallyzmin = math.floor(zmin / zstep) * zstep
tallyzmax = math.ceil(zmax / zstep) * zstep 
xdim = int(abs(tallyxmax - tallyxmin) / xstep)
ydim = int(abs(tallyymax - tallyymin) / ystep)
zdim = int(abs(tallyzmax - tallyzmin) / zstep)

regularmesh = openmc.RegularMesh(name = "Regular Mesh")
regularmesh.dimension = (xdim, ydim, zdim)
regularmesh.lower_left = (tallyxmin, tallyymin, tallyzmin)
#regularmesh.upper_right = (tallyxmax, tallyymax, tallyzmax)
regularmesh.width = (xstep, ystep, zstep)
print(regularmesh)

xstep = 2 * hw
ystep = 2 * hw
zstep = 100
tallyzmin = math.floor(zmin / zstep) * zstep
tallyzmax = math.ceil(zmax / zstep) * zstep 
xdim = 1
ydim = 1
zdim = int(abs(tallyzmax - tallyzmin) / zstep)

zmesh = openmc.RegularMesh(name = "z-direction Mesh")
zmesh.dimension = (xdim, ydim, zdim)
zmesh.lower_left = (xmin, ymin, tallyzmin)
#zmesh.upper_right = (tallyxmax, tallyymax, tallyzmax)
zmesh.width = (xstep, ystep, zstep)

energyfilter = openmc.EnergyFilter(np.logspace(-3, 8, num=100))
meshfilter = openmc.MeshFilter(regularmesh)
zmeshfilter = openmc.MeshFilter(zmesh)
meshsurfacefilter = openmc.MeshSurfaceFilter(regularmesh)

tally = openmc.Tally(name = "flux, regular mesh")
tally.filters = [meshfilter]
tally.scores = ["flux"]
tallies.append(tally)

tally = openmc.Tally(name = "flux, z-direction mesh")
tally.filters = [zmeshfilter, energyfilter]
tally.scores = ["flux"]
tallies.append(tally)

tally = openmc.Tally(name = "current, regular mesh")
tally.filters = [meshsurfacefilter]
tally.scores = ["current"]
tallies.append(tally)

for p in [concretepath, waterpath, soilpath,
          concrete2100path, water2100path, soil2100path,
          concretefetterpath, waterfetterpath, soilfetterpath]:
    print("Tallies", p)
    tallies.export_to_xml(os.path.join(p, "tallies.xml"))
