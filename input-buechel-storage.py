import openmc
import json
import os
import matplotlib.pyplot as plt
import math
import numpy as np
import sys
import argparse

import helper

cspath = "/openmc/openmc-data/v0.12/lib80x_hdf5/cross_sections.xml"

basepath = "bust-20210714-fetter"
fettersource = False
particles = 1
batches = 10
plot = False
weaponage = 0 # years

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--fetter", action="store_true",
                    help="add fetter model as source")
parser.add_argument("-p", "--plot", action="store_true",
                    help="plot some cuts through the geometry")
parser.add_argument("-b", "--basepath", 
                    help="set basepath for calculation")
args = parser.parse_args()
if args.fetter:
    fettersource = True
if args.plot:
    plot = True
if not args.basepath is None:
    basepath = args.basepath

if not os.path.exists(basepath):
    os.mkdir(basepath)

# dump basic settings to json file
with open(os.path.join(basepath, "calculation.json"), 'w') as f:
    json.dump({'particles': particles,
               'basepath': basepath,
               'cspath': cspath,
               'batches': batches,
               'fettersource': fettersource,
               'weaponage': weaponage},
              f)
    f.close()

###############################################################################
# Materials
###############################################################################

steelMat = openmc.Material(name = "Steel (for doors)")
steelMat.set_density("g/cm3", 7.82) # PNNL Compendium page 281
steelMat.add_element("C", 0.022831, "ao")
steelMat.add_element("Fe", 0.977169, "ao")

concreteRegularMat = openmc.Material(name = "Regular Concrete")
concreteRegularMat.set_density("g/cm3", 2.3) # PNNL Compendium page 112
concreteRegularMat.add_element("H", 0.16038, "ao")
concreteRegularMat.add_element("O", 0.563183, "ao")
concreteRegularMat.add_element("Na", 0.021365, "ao")
concreteRegularMat.add_element("Al", 0.021343, "ao")
concreteRegularMat.add_element("Si", 0.203231, "ao")
concreteRegularMat.add_element("Ca", 0.018595, "ao")
concreteRegularMat.add_element("Fe", 0.004246, "ao")

# 300kg carbon steel / m3
# https://www.yourspreadsheets.co.uk/reinforcement-estimates.html
steelweightperm3 = 300
steelvolumefraction = 300 / (7.82 * 1e4)
concretevolumefraction = 1 - steelvolumefraction
wallMat = openmc.Material.mix_materials([concreteRegularMat, steelMat],
                                        [concretevolumefraction, steelvolumefraction], 'vo')
wallMat.name = "Wall Material"

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

#*******************************************************************************
# fetter model

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
ages = [weaponage]
puvec.createagedvector(ages)

pudf = puvec.pudf
puagedf = puvec.puagedf

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

materiallist = [steelMat, concreteRegularMat, wallMat, airMat, soilMat]
if fettersource:
    materiallist = materiallist + [wpuMat[weaponage], berMat, tunMat, hmxMat, aluMat]

materials = openmc.Materials(materiallist)
materials.cross_sections = cspath
materials.export_to_xml(os.path.join(basepath, "materials.xml"))

###############################################################################
# Geometry
###############################################################################

outermarginx = 2000
outermarginy = 2000
outermarginz = 500

vaulthorizontalpixel = 359.7 / 564.7 # from drawing - warhead 11.8 feet / 359.7cm
weapondiameter = 33
vaultverticalpixel = weapondiameter / 62.1 # from drawing - warhead 12 inches / 33 cm
vaultt1 = 60.1 * vaultverticalpixel
vaultt2 = 23.8 * vaultverticalpixel
vaultd1 = 411 * vaultverticalpixel
vaultd2 = (vaultd1 - vaultt1 - vaultt2) / 2
vaultd3 = 2.5 * weapondiameter + 2 * vaultt2
vaultd4 = 655.9 * vaulthorizontalpixel
vaultd5 = vaultd4 + 2 * vaultt2 + 2 * vaultt1

paspixel = 2770 / 304.1
pasw1 = (164.7 / 2) * paspixel
pasw2 = (304.1 / 2) * paspixel
pash1 = 12.1 * paspixel
pash2 = 33.2 * paspixel
pash3 = 55.7 * paspixel
pasr1 = (pasw1 ** 2 + pash2 ** 2) / (2 *  pash2)
pasr2 = pasr1 - pash1
pash4 = pasr1 - (pash1 + pash3 + pash2)
pasw3 = pasw1 - pash1 - vaultd5

basepixel = 30000 / 981.4
pasd1 = 100.7 * basepixel
pasd2 = 300 # estimate
past1 = 10 # steel door thickness
past2 = pash1
pasd3 = pasd1 - (pasd2 + past1 + past2)
pasd4 = pasd2 + past2 + 500 # from Kristensen 2005 figure (15 feet)

surfaces = {}

# x surfaces for pas
xmin = -pasw2 - outermarginx
xmax =  pasw2 + outermarginx

surfacex = { 1: xmin, # never set surface id to 0 - then c function can't find cell
             10: -pasw2,
             20: -pasw1,
             30: 0,
             40: pasw1,
             50: pasw2,
             60: xmax,
             }

for sid, val in surfacex.items():
    surfaces[sid] = openmc.XPlane(val, surface_id = sid)

surfaces[1].boundary_type = 'vacuum'
surfaces[60].boundary_type = 'vacuum'

# z surfaces + cylinders
zmin = 0
zmax = pasr1 + outermarginz

surfaces[2000] = openmc.ZPlane(zmin, surface_id = 2000)
surfaces[2010] = openmc.ZPlane(pash4, surface_id = 2010)
surfaces[2020] = openmc.ZPlane(pash4 + pash1, surface_id = 2020)

surfaces[2030] = openmc.YCylinder(r=pasr2, surface_id = 2030)
surfaces[2040] = openmc.YCylinder(r=pasr1, surface_id = 2040)

surfaces[2050] = openmc.ZPlane(zmax, surface_id = 2050)

# left inclined plane
a, b, c, d = helper.planeparameter3points([-pasw2, 0, pash4 + pash1],
                                   [-pasw2, pasd2, pash4 + pash1],
                                   [-pasw1, 0, pash4 + pash1 + pash3])
surfaces[2060] = openmc.Plane(a, b, c, d, surface_id = 2060)

# right inclined plane
a, b, c, d = helper.planeparameter3points([pasw2, 0, pash4 + pash1],
                                   [pasw1, 0, pash4 + pash1 + pash3],
                                   [pasw2, pasd2, pash4 + pash1])
surfaces[2070] = openmc.Plane(a, b, c, d, surface_id = 2070)

surfaces[2000].boundary_type = 'vacuum'
surfaces[2050].boundary_type = 'vacuum'
    
# y surfaces pas

ymin = -outermarginy
ymax = pasd1 + outermarginy

surfaces[1000] = openmc.YPlane(ymin, surface_id = 1000)

surfaces[1010] = openmc.YPlane(0, surface_id = 1010)
surfaces[1020] = openmc.YPlane(pasd2, surface_id = 1020)
surfaces[1030] = openmc.YPlane(pasd2 + past1, surface_id = 1030)
surfaces[1040] = openmc.YPlane(pasd1 - pash1, surface_id = 1040)
surfaces[1050] = openmc.YPlane(pasd1, surface_id = 1050)

surfaces[1060] = openmc.YPlane(ymax, surface_id = 1060)

surfaces[1000].boundary_type = 'vacuum'
surfaces[1060].boundary_type = 'vacuum'


# x surfaces vault
xoffset = pasw3
surfaces[3000] = openmc.XPlane(xoffset, surface_id = 3000)
surfaces[3010] = openmc.XPlane(xoffset + vaultt1, surface_id = 3010)
surfaces[3020] = openmc.XPlane(xoffset + vaultt1 + vaultt2, surface_id = 3020)
surfaces[3030] = openmc.XPlane(xoffset + vaultt1 + 2 * vaultt2, surface_id = 3030)
surfaces[3040] = openmc.XPlane(xoffset + vaultt1 + vaultd4, surface_id = 3040)
surfaces[3050] = openmc.XPlane(xoffset + vaultt1 + vaultt2 + vaultd4, surface_id = 3050)
surfaces[3060] = openmc.XPlane(xoffset + vaultt1 + vaultt2 + vaultd4 + vaultt2, surface_id = 3060)
surfaces[3070] = openmc.XPlane(xoffset + vaultd5, surface_id = 3070)

# y surfaces vault
yoffset = pasd4
surfaces[3100] = openmc.YPlane(yoffset, surface_id = 3100)
surfaces[3110] = openmc.YPlane(yoffset + vaultt1, surface_id = 3110)
surfaces[3120] = openmc.YPlane(yoffset + vaultt1 + vaultt2, surface_id = 3120)
surfaces[3130] = openmc.YPlane(yoffset + vaultt1 + 2 * vaultt2, surface_id = 3130)
surfaces[3140] = openmc.YPlane(yoffset + vaultt1 + vaultd3 - 2 * vaultt2, surface_id = 3140)
surfaces[3150] = openmc.YPlane(yoffset + vaultt1 + vaultd3 - vaultt2, surface_id = 3150)
surfaces[3160] = openmc.YPlane(yoffset + vaultt1 + vaultd3, surface_id = 3160)
surfaces[3170] = openmc.YPlane(yoffset + vaultt1 + vaultd3 + vaultt1, surface_id = 3170)

# z surfaces vault
zoffset = pash4 + pash1
surfaces[3250] = openmc.ZPlane(zoffset, surface_id = 3250)
surfaces[3240] = openmc.ZPlane(zoffset - vaultt1, surface_id = 3240)
surfaces[3230] = openmc.ZPlane(zoffset - vaultt1 - vaultd2, surface_id = 3230)
surfaces[3220] = openmc.ZPlane(zoffset - vaultt1 - vaultd2 - vaultt2, surface_id = 3220)
surfaces[3210] = openmc.ZPlane(zoffset - vaultd1, surface_id = 3210)
surfaces[3200] = openmc.ZPlane(zoffset - vaultd1 - vaultt1, surface_id = 3200)

vaultregion = +surfaces[3000] & -surfaces[3070] & +surfaces[3100] & -surfaces[3170] & +surfaces[3200] & -surfaces[3250]

pascells = []
soilcell = openmc.Cell()
soilcell.region = \
    +surfaces[10] & -surfaces[50] &\
    +surfaces[1010] & -surfaces[1050] &\
    +surfaces[2000] & -surfaces[2010] &\
   ~vaultregion
soilcell.fill = soilMat
soilcell.name = "Soil Layer"
pascells.append(soilcell)

concretefloor = openmc.Cell()
concretefloor.region = \
    +surfaces[10] & -surfaces[50] &\
    +surfaces[1010] & -surfaces[1050] &\
    +surfaces[2010] & -surfaces[2020] &\
   ~vaultregion
concretefloor.fill = concreteRegularMat # no steel here
concretefloor.name = "Concrete Floor Layer"
pascells.append(concretefloor)

roofcell = openmc.Cell()
roofcell.region = \
    +surfaces[20] & -surfaces[40] &\
    +surfaces[1010] & -surfaces[1050] &\
    +surfaces[2030] & -surfaces[2040]
roofcell.fill = wallMat
roofcell.name = "PAS Roof"
pascells.append(roofcell)

leftwall = openmc.Cell()
leftwall.region = \
    +surfaces[2020] & -surfaces[2060] & - surfaces[20] &\
    +surfaces[1010] & -surfaces[1050]
leftwall.fill = wallMat
leftwall.name = "PAS Left Wall"
pascells.append(leftwall)

rightwall = openmc.Cell()
rightwall.region = \
    +surfaces[2020] & -surfaces[2070] & + surfaces[40] &\
    +surfaces[1010] & -surfaces[1050]
rightwall.fill = wallMat
rightwall.name = "PAS Right Wall"
pascells.append(rightwall)

backwall = openmc.Cell()
backwall.region = \
    +surfaces[20] & -surfaces[40] &\
    +surfaces[1040] & -surfaces[1050] &\
    +surfaces[2020] & -surfaces[2030]
backwall.fill = wallMat
backwall.name = "PAS Back Wall"
pascells.append(backwall)

pasdoor = openmc.Cell()
pasdoor.region =\
    +surfaces[20] & -surfaces[40] &\
    +surfaces[1020] & -surfaces[1030] &\
    +surfaces[2020] & -surfaces[2030]
pasdoor.fill = steelMat
pasdoor.name = "PAS Main Door"
pascells.append(pasdoor)

pasinnerair = openmc.Cell()
pasinnerair.region =\
    +surfaces[20] & -surfaces[40] &\
    +surfaces[1010] & -surfaces[1040] &\
    +surfaces[2020] & -surfaces[2030] &\
    ~pasdoor.region
pasinnerair.fill = airMat
pasinnerair.name = "PAS Air"
pascells.append(pasinnerair)

pasouterair = openmc.Cell()
pasouterair.region =\
    +surfaces[10] & -surfaces[50] &\
    +surfaces[1010] & -surfaces[1050] &\
    +surfaces[2020] & +surfaces[2040] & -surfaces[2050] &\
    ~leftwall.region & ~rightwall.region
pasouterair.fill = airMat
pasouterair.name = "PAS Outer Air"
pascells.append(pasouterair)

vaultcells = []
vaultairregion = +surfaces[3010] & -surfaces[3060] & +surfaces[3110] & -surfaces[3160] & +surfaces[3210] & -surfaces[3240]

vaultconcrete = openmc.Cell()
vaultconcrete.region =\
    +surfaces[3000] & -surfaces[3070] &\
    +surfaces[3100] & -surfaces[3170] &\
    +surfaces[3200] & -surfaces[3250] &\
    ~vaultairregion
vaultconcrete.fill = wallMat
vaultconcrete.name = "Vault Concrete Outer Layer"
vaultcells.append(vaultconcrete)

vaultpoles = {i: openmc.Cell() for i in range(1, 5)}
vaultpoles[1].region =\
    +surfaces[3020] & -surfaces[3030] &\
    +surfaces[3120] & -surfaces[3130] &\
    +surfaces[3210] & -surfaces[3240]
vaultpoles[2].region =\
    +surfaces[3020] & -surfaces[3030] &\
    +surfaces[3140] & -surfaces[3150] &\
    +surfaces[3210] & -surfaces[3240]
vaultpoles[3].region =\
    +surfaces[3040] & -surfaces[3050] &\
    +surfaces[3120] & -surfaces[3130] &\
    +surfaces[3210] & -surfaces[3240]
vaultpoles[4].region =\
    +surfaces[3040] & -surfaces[3050] &\
    +surfaces[3140] & -surfaces[3150] &\
    +surfaces[3210] & -surfaces[3240]

for i in vaultpoles:
    vaultpoles[i].fill = steelMat
    vaultpoles[i].name = "Vault Pole"
    vaultcells.append(vaultpoles[i])

vaultshelf = openmc.Cell()
vaultshelf.region =\
    +surfaces[3030] & -surfaces[3040] &\
    +surfaces[3120] & -surfaces[3150] &\
    +surfaces[3220] & -surfaces[3230]
vaultshelf.fill = steelMat
vaultshelf.name = "Vault Shelf"
vaultcells.append(vaultshelf)

vaultair = openmc.Cell()
vaultair.region =\
    vaultairregion & ~vaultshelf.region &\
    ~vaultpoles[1].region & ~vaultpoles[2].region &\
    ~vaultpoles[3].region & ~vaultpoles[4].region
vaultair.fill = airMat
vaultair.name = "Vault Air"
#vaultcells.append(vaultair)

margincells = []
pasregion = +surfaces[10] & -surfaces[50] & +surfaces[1010] & -surfaces[1050]

marginsoil = openmc.Cell()
marginsoil.region = \
    +surfaces[1] & -surfaces[60] &\
    +surfaces[1000] & -surfaces[1060] &\
    +surfaces[2000] & -surfaces[2020] &\
    ~pasregion
marginsoil.fill = soilMat
marginsoil.name = "Margin - Soil"
margincells.append(marginsoil)

marginair = {i: openmc.Cell() for i in range(1, 5)}
marginair[1].region = \
    +surfaces[1] & -surfaces[60] &\
    +surfaces[1000] & -surfaces[1010] &\
    +surfaces[2020] & -surfaces[2050]
marginair[2].region = \
    +surfaces[1] & -surfaces[60] &\
    +surfaces[1050] & -surfaces[1060] &\
    +surfaces[2020] & -surfaces[2050]
marginair[3].region = \
    +surfaces[1] & -surfaces[10] &\
    +surfaces[1010] & -surfaces[1050] &\
    +surfaces[2020] & -surfaces[2050]
marginair[4].region = \
    +surfaces[50] & -surfaces[60] &\
    +surfaces[1010] & -surfaces[1050] &\
    +surfaces[2020] & -surfaces[2050]

for i in [1, 2, 3, 4]:
    marginair[i].fill = airMat
    marginair[i].name = "Margin - Air"
    margincells.append(marginair[i])


#*******************************************************************************
# fetter model

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

airCell = openmc.Cell()
airCell.region = +aluSurface
airCell.name = "Air around weapon"

pitCell.fill = wpuMat[weaponage] 
refCell.fill = berMat
tamCell.fill = tunMat
expCell.fill = hmxMat
aluCell.fill = aluMat
airCell.fill = airMat

#*******************************************************************************
# cells & universes + source location

if fettersource:
    # Source in lower vault compartment, at center line of one weapon
    sourcex = surfaces[3000].x0 + vaultd5 / 4
    # fetter model in center of vault
    sourcey = surfaces[3110].y0 + vaultd3 / 4
    sourcez = surfaces[3210].z0 + vaultd2 / 2
    # make space in vault
    weaponsurface = openmc.Sphere(r = aluOR + 1, x0 = sourcex, y0 = sourcey, z0 = sourcez) # slightly bigger to avoid particles getting 'trapped' between universes, reaching MAX_EVENTS
    vaultair.region &= +weaponsurface
else:
    # Source in lower vault compartment, at center line of one weapon
    sourcex = surfaces[3000].x0 + vaultd5 / 4
    sourcey = surfaces[3110].y0 + vaultt2 + weapondiameter / 2
    sourcez = surfaces[3210].z0 + vaultd2 / 2

vaultcells.append(vaultair)
cells = pascells + vaultcells + margincells
root = openmc.Universe(cells = cells)

if fettersource:
    fettercells = [cenCell, pitCell, refCell, tamCell, expCell, aluCell, airCell]
    fetteruniverse = openmc.Universe(cells = fettercells)

    weaponcell = openmc.Cell(fill = fetteruniverse, region = -weaponsurface)
    weaponcell.translation = (sourcex, sourcey, sourcez)
    root.add_cell(weaponcell)

geometry = openmc.Geometry(root)
geometry.export_to_xml(os.path.join(basepath, "geometry.xml"))

#*******************************************************************************
# Plot
if plot:
    print("Plotting geometry - can take a few seconds")
    xfactor = 1.2
    yfactor = 1.2
    zfactor = 1.2

    plt.figure(figsize=(8, 8))
    root.plot(origin = (sourcex, sourcey, sourcez),
              basis=('xz'),
              width=(300, 300),
              pixels=(600, 600),
              seed = 1)
    plt.title("Source View (xz)")
    plt.savefig(os.path.join(basepath, "source-view-xz.png"))
    plt.close()

    plt.figure(figsize=(8, 8))
    root.plot(origin = (sourcex, sourcey, sourcez),
              basis=('xy'),
              width=(300, 300),
              pixels=(600, 600),
              seed = 1)
    plt.title("Source View (xy)")
    plt.savefig(os.path.join(basepath, "source-view-xy.png"))
    plt.close()
    
    plt.figure(figsize=(8, 8))
    root.plot(origin = (0, 1, (zmax - zmin) / 2), basis=('xz'), width=((xmax - xmin) * xfactor, (zmax - zmin) * zfactor), seed = 1)
    plt.title("Front View")
    plt.savefig(os.path.join(basepath, "front-view.png"))
    plt.close()

    plt.figure(figsize=(8, 8))
    root.plot(origin = (0, pasd4 + vaultt1 + vaultd3 / 2, (zmax - zmin) / 2), basis=('xz'), width=((xmax - xmin) * xfactor, (zmax - zmin) * zfactor), seed = 1)
    plt.title("Front View, Vault position")
    plt.savefig(os.path.join(basepath, "front-view-vault.png"))
    plt.close()

    plt.figure(figsize=(8, 8))
    root.plot(origin = (0, ymin + (ymax - ymin) / 2, pash4 + pash1 + 1),
              basis=('xy'),
              width=((xmax - xmin) * xfactor, (ymax - ymin) * yfactor))
    plt.title("Top View, Floor Level")
    plt.savefig(os.path.join(basepath, "top-view-floor-level.png"))
    plt.close()

    for z in range(math.ceil(zmin / 100), math.ceil(zmax / 100)):
        plt.figure(figsize=(8, 8))
        root.plot(origin = (0, ymin + (ymax - ymin) / 2, z * 100),
                  basis=('xy'),
                  width=((xmax - xmin) * xfactor, (ymax - ymin) * yfactor),
                  seed = 1)
        plt.title("Top View, Floor Level, z = {:03d}m".format(z))
        plt.savefig(os.path.join(basepath, "top-view-floor-level-{:03d}.png".format(z)))
        plt.close()

    plt.figure(figsize=(8, 8))
    root.plot(origin = (pasw3 + vaultd5 / 2, pasd4 + vaultt1 + vaultt2 * 1.5, pash4 + pash1 - (vaultd1 + vaultt1) / 2), basis=('xz'), width=(vaultd5 * xfactor, (vaultd1 + vaultt1) * zfactor), seed = 1)
    plt.title("Vault Front View")
    plt.savefig(os.path.join(basepath, "vault-front-view.png"))
    plt.close()


###############################################################################
# Source & Settings
###############################################################################

if fettersource:
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
else:
    source = openmc.Source(space = openmc.stats.Point((sourcex, sourcey, sourcez)),
                           energy = openmc.stats.Discrete([2e6], [1]),
                           particle = 'neutron')

settings = openmc.Settings()
settings.run_mode = 'fixed source'
settings.inactive = 0
settings.batches = batches
settings.particles = particles
settings.source = source

settings.export_to_xml(os.path.join(basepath, "settings.xml"))

###############################################################################
# Tallies
###############################################################################

tallies = openmc.Tallies()

xwidth = xmax - xmin
ywidth = ymax - ymin
zwidth = zmax - zmin

xstep = 100
tallyxmin = math.floor(xmin / xstep) * xstep
tallyxmax = math.ceil(xmax / xstep) * xstep + 1
ystep = 100
tallyymin = math.floor(ymin / ystep) * ystep
tallyymax = math.ceil(ymax / ystep) * ystep + 1
zstep = 100
tallyzmin = math.floor(zmin / zstep) * zstep
tallyzmax = math.ceil(zmax / zstep) * zstep + 1

mainmesh = openmc.RectilinearMesh(name = "Main Mesh")
mainmesh.x_grid = range(tallyxmin, tallyxmax, xstep)
mainmesh.y_grid = range(tallyymin, tallyymax, ystep)
mainmesh.z_grid = range(tallyzmin, tallyzmax, zstep)

meshfilter = openmc.MeshFilter(mainmesh)
tally = openmc.Tally(name = "flux")
tally.filters = [meshfilter]
tally.scores = ["flux"]
tallies.append(tally)

tallies.export_to_xml(os.path.join(basepath, "tallies.xml"))
