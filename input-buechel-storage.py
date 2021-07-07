import openmc
import json
import os
import matplotlib.pyplot as plt
import math
import numpy as np
import sys

basepath = "bust-20210707"
particles = 1

cspath = "/openmc/openmc-data/v0.12/lib80x_hdf5/cross_sections.xml"

if not os.path.exists(basepath):
    os.mkdir(basepath)

# dump basic settings to json file
with open(os.path.join(basepath, "calculation.json"), 'w') as f:
    json.dump({'particles': particles,
               'basepath': basepath,
               'cspath': cspath},
              f)
    f.close()

def planeparameter3points(p1, p2, p3):
    if len(p1) != 3 and len(p2) != 3 and len(p3) != 3:
        print("All points need 3 dimensions")
        raise RuntimeError()
    p1 = np.array(p1)
    p2 = np.array(p2)
    p3 = np.array(p3)

    v1 = p3 - p1
    v2 = p2 - p1

    cp = np.cross(v1, v2)
    a, b, c = cp

    d = np.dot(cp, p3)

    return (a, b, c, d)

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

materiallist = [steelMat, concreteRegularMat, wallMat, airMat, soilMat]
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
a, b, c, d = planeparameter3points([-pasw2, 0, pash4 + pash1],
                                   [-pasw2, pasd2, pash4 + pash1],
                                   [-pasw1, 0, pash4 + pash1 + pash3])
surfaces[2060] = openmc.Plane(a, b, c, d, surface_id = 2060)

# right inclined plane
a, b, c, d = planeparameter3points([pasw2, 0, pash4 + pash1],
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
vaultcells.append(vaultair)

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

cells = pascells + vaultcells + margincells
#cells = margincells

root = openmc.Universe(cells = cells)
geometry = openmc.Geometry(root)
geometry.export_to_xml(os.path.join(basepath, "geometry.xml"))

# Plot
if len(sys.argv) > 1:
    print("Plotting geometry - can take a few seconds")
    xfactor = 1.2
    yfactor = 1.2
    zfactor = 1.2
    plt.figure(figsize=(8, 8))
    root.plot(origin = (0, 1, (zmax - zmin) / 2), basis=('xz'), width=((xmax - xmin) * xfactor, (zmax - zmin) * zfactor))
    plt.title("Front View")
    plt.savefig(os.path.join(basepath, "front-view.png"))

    plt.figure(figsize=(8, 8))
    root.plot(origin = (0, pasd4 + vaultt1 + vaultd3 / 2, (zmax - zmin) / 2), basis=('xz'), width=((xmax - xmin) * xfactor, (zmax - zmin) * zfactor))
    plt.title("Front View, Vault position")
    plt.savefig(os.path.join(basepath, "front-view-vault.png"))

    plt.figure(figsize=(8, 8))
    root.plot(origin = (0, ymin + (ymax - ymin) / 2, pash4 + pash1 + 1),
              basis=('xy'),
              width=((xmax - xmin) * xfactor, (ymax - ymin) * yfactor))
    plt.title("Top View, Floor Level")
    plt.savefig(os.path.join(basepath, "top-view-floor-level.png"))

    for z in range(math.ceil(zmin / 100), math.ceil(zmax / 100)):
        plt.figure(figsize=(8, 8))
        root.plot(origin = (0, ymin + (ymax - ymin) / 2, z * 100),
                  basis=('xy'),
                  width=((xmax - xmin) * xfactor, (ymax - ymin) * yfactor))
        plt.title("Top View, Floor Level, z = {:03d}m".format(z))
        plt.savefig(os.path.join(basepath, "top-view-floor-level-{:03d}.png".format(z)))

    plt.figure(figsize=(8, 8))
    root.plot(origin = (pasw3 + vaultd5 / 2, pasd4 + vaultt1 + vaultt2 * 1.5, pash4 + pash1 - (vaultd1 + vaultt1) / 2), basis=('xz'), width=(vaultd5 * xfactor, (vaultd1 + vaultt1) * zfactor))
    plt.title("Vault Front View")
    plt.savefig(os.path.join(basepath, "vault-front-view.png"))


###############################################################################
# Source & Settings
###############################################################################

# Source in lower vault compartment, at center line of one weapon
sourcex = surfaces[3000].x0 + vaultd5 / 4
sourcey = surfaces[3110].y0 + vaultt2 + weapondiameter / 2
sourcez = surfaces[3210].z0 + vaultd2 / 2

source = openmc.Source(space = openmc.stats.Point((sourcex, sourcey, sourcez)),
                       energy = openmc.stats.Discrete([2e6], [1]),
                       particle = 'neutron')


settings = openmc.Settings()
settings.run_mode = 'fixed source'
settings.inactive = 0
settings.batches = 10
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
xmax = math.ceil(xwidth / 100) * 100 + 1
ystep = 100
ymax = math.ceil(ywidth / 100) * 100 + 1
zstep = 100
zmax = math.ceil(zwidth / 100) * 100 + 1

mainmesh = openmc.RectilinearMesh(name = "Main Mesh")
mainmesh.x_grid = range(0, xmax, xstep)
mainmesh.y_grid = range(0, ymax, ystep)
mainmesh.z_grid = range(0, zmax, zstep)

meshfilter = openmc.MeshFilter(mainmesh)
tally = openmc.Tally(name = "flux")
tally.filters = [meshfilter]
tally.scores = ["flux"]
tallies.append(tally)

tallies.export_to_xml(os.path.join(basepath, "tallies.xml"))
