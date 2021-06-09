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

basepath = "20210605-100k-evonly"
particles = 100000

cspath = "/openmc/openmc-data/v0.12/lib80x_hdf5/cross_sections.xml"
ages = [0, 5, 10, 15, 20, 25, 30, 35, 40]

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
        os.mkdir(os.path.join(evpath, "full-{:02d}".format(age)))
        os.mkdir(os.path.join(fspath, "full-{:02d}".format(age)))

# dump basic settings to json file
with open(os.path.join(basepath, "calculation.json"), 'w') as f:
    json.dump({'ages': ages,
               'geometries': geometries,
               'particles': particles,
               'basepath': basepath,
               'cspath': cspath},
              f)
    f.close()
    
###############################################################################
# Materials
###############################################################################

# need pit radius already here
pitOR = 5
pitThickness = 0.75
pitIR = pitOR - pitThickness


pudict = {'name': ['Pu238', 'Pu239', 'Pu240', 'Pu241', 'Pu242', 'Am241'],
          'wo': [0.00005, 0.933, 0.060, 0.0044, 0.00015, 0], # Fetter et al. 1990
          'sfneutrons': [2630, 0.0152, 1031, 0.001723, 1722, 1.47], # PhD Thesis Moritz
          'watt-a': [1.17948e-6, 1.12963e-6, 1.25797e-6, 1.18698e-6, 1.22078e-6, 1.07179e-6], # Verbeke et al. 2010
          'watt-b': [4.16933e-6, 3.80269e-6, 4.68927e-6, 4.15150e-6, 4.36668e-6, 3.46195e-6]} # Verbeke et al. 2010
pudict['mass'] = [openmc.data.atomic_mass(iso) for iso in pudict['name']]
pudf = pd.DataFrame(pudict)
pudf.set_index('name', inplace = True)

puagedf = pd.DataFrame(columns = ['name', 'age', 'wo', 'sfneutrons', 'watt-a', 'watt-b'])

# all times in years
lbda = np.log(2)/14.329 # T_1/2 from PhD Thesis Moritz
totmol = sum(pudf['wo'] / pudf['mass']) + 0.002 / openmc.data.atomic_weight('O')

Npu0 = pudf.loc['Pu241', 'wo'] / openmc.data.atomic_mass('Pu241') / totmol 
for age in ages:
    agelist = [age] * 6
    Npu = Npu0 * np.exp(-lbda * age)
    Nam = Npu0 * (1 - np.exp(-lbda * age))
    puadict = copy.deepcopy(pudict)
    puadict['wo'][3] = Npu * totmol * openmc.data.atomic_mass('Pu241')
    puadict['wo'][5] = Nam * totmol * openmc.data.atomic_mass('Am241')
    puadict['age'] = agelist

    puagedf = pd.concat([puagedf, pd.DataFrame(puadict)])
    
puagedf.set_index(['name', 'age'], inplace=True)

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
    materials.export_to_xml(os.path.join(evpath, "full-{:02d}".format(age), "materials.xml"))
    materials.export_to_xml(os.path.join(fspath, "full-{:02d}".format(age), "materials.xml"))


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
    geometry.export_to_xml(os.path.join(evpath, "full-{:02d}".format(age), "geometry.xml"))
    geometry.export_to_xml(os.path.join(fspath, "full-{:02d}".format(age), "geometry.xml"))

plt.figure(figsize=(5, 5))
root.plot(width=(47, 47))
plt.title("Geometry of Fetter model in OpenMC")
plt.savefig(os.path.join(basepath, "fetter_geometry.png"))


###############################################################################
# Tallies
###############################################################################

tallies = openmc.Tallies()

# Flux
tally = openmc.Tally(name='flux')
tally.filters = [openmc.CellFilter([pitCell, refCell, tamCell, expCell, aluCell, outCell]),
                 openmc.EnergyFilter(np.logspace(-3, 8, num=100))]
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

for geo in geometries:
    tallies.export_to_xml(os.path.join(fspath, geo, "tallies.xml"))
    tallies.export_to_xml(os.path.join(evpath, geo, "tallies.xml"))
for age in ages:
    tallies.export_to_xml(os.path.join(fspath, "full-{:02d}".format(age)))
    tallies.export_to_xml(os.path.join(evpath, "full-{:02d}".format(age)))
                 
###############################################################################
# Source & Settings
###############################################################################


bounds = [-pitOR, -pitOR, -pitOR, pitOR, pitOR, pitOR]
uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable=True)
src = openmc.Source(space = uniform_dist, particle = 'neutron')

settings = openmc.Settings()
settings.run_mode = 'eigenvalue'
settings.batches = 120
settings.inactive = 20
settings.particles = particles
settings.source = src

for geo in geometries:
    settings.export_to_xml(os.path.join(fspath, geo, "settings.xml"))
    settings.export_to_xml(os.path.join(evpath, geo, "settings.xml"))
for age in ages:
    settings.export_to_xml(os.path.join(fspath, "full-{:02d}".format(age)))
    settings.export_to_xml(os.path.join(evpath, "full-{:02d}".format(age)))


#openmc.data.atomic_mass('Fe54')
