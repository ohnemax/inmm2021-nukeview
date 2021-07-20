import json
import os
import matplotlib.pyplot as plt
import math
import numpy as np
import sys
import argparse

import openmc
import helper

cspath = "/openmc/openmc-data/v0.12/lib80x_hdf5/cross_sections.xml"

basepath = "jovanovic-erickson-example"
particles = 10000
batches = 10
plot = False
survival = False
time = 0

parser = argparse.ArgumentParser()
parser.add_argument("-p", "--plot", action="store_true",
                    help="plot some cuts through the geometry")
parser.add_argument("-s", "--survival", action="store_true",
                    help="turn on survival biasing")
parser.add_argument("-b", "--basepath", 
                    help="set basepath for calculation")
parser.add_argument("-c", "--batches", type = int,
                    help="set number of batches")
parser.add_argument("-n", "--particles", type = int,
                    help="set number of particles")
parser.add_argument("-t", "--time", type = float,
                    help="set a measurement time (and with that the particle number)")

args = parser.parse_args()
if args.plot:
    plot = True
if args.survival:
    survival = True
if not args.basepath is None:
    basepath = args.basepath
if not args.batches is None:
    batches = args.batches
if not args.particles is None:
    particles = args.particles
if not args.time is None:
    time = args.time

if not os.path.exists(basepath):
    os.mkdir(basepath)

# dump basic settings to json file
with open(os.path.join(basepath, "calculation.json"), 'w') as f:
    json.dump({'particles': particles,
               'basepath': basepath,
               'cspath': cspath,
               'batches': batches},
              f)
    f.close()

###############################################################################
# Materials
###############################################################################

puvec = helper.puvector()
puvec.setwo('Pu238', 0.00005)
puvec.setwo('Pu239', 0.993)
puvec.setwo('Pu240', 0.060)
puvec.setwo('Pu241', 0.0044)
puvec.setwo('Pu242', 0.00015)
puvec.setwo('Am241', 0) # from fetter-1990
puvec.createagedvector([0])
pudf = puvec.pudf

wpuMat = openmc.Material(name = "Weapon-grade Plutonium")
wpuMat.set_density("g/cm3", 19.84) # PNNL Compendium, page 233
wpuMat.add_nuclide("Pu238", pudf.loc['Pu238', 'wo'], "wo")
wpuMat.add_nuclide("Pu239", pudf.loc['Pu239', 'wo'], "wo")
wpuMat.add_nuclide("Pu240", pudf.loc['Pu240', 'wo'], "wo") #should be wo all along
wpuMat.add_nuclide("Pu241", pudf.loc['Pu241', 'wo'], "wo") # need to consider for decay!
wpuMat.add_nuclide("Pu242", pudf.loc['Pu242', 'wo'], "wo")
wpuMat.add_nuclide("Am241", pudf.loc['Am241', 'wo'], "wo")
wpuMat.add_element("O", 0.002, "wo")

watMat = openmc.Material(name = "Water")
watMat.set_density("g/cm3", 0.998207) #PNNL Compendium page 338
watMat.add_element("H", 2, "ao")
watMat.add_element("O", 1, "ao")
watMat.add_s_alpha_beta('c_H_in_H2O')

materiallist = [wpuMat, watMat]

materials = openmc.Materials(materiallist)
materials.cross_sections = cspath
materials.export_to_xml(os.path.join(basepath, "materials.xml"))

###############################################################################
# Geometry
###############################################################################

wpumass = 5000

wpuvolume = wpumass / wpuMat.density
wpuradius = (3 / 4 / np.pi * wpuvolume) ** (1 / 3)

wpuSurface = openmc.Sphere(r = wpuradius)

hw = 100
watTop = openmc.ZPlane(hw)
watBottom = openmc.ZPlane(-hw)
watLeft = openmc.XPlane(-hw)
watRight = openmc.XPlane(hw)
watFront = openmc.YPlane(-hw)
watBack = openmc.YPlane(hw)

vacTop = openmc.ZPlane(2 * hw)
vacBottom = openmc.ZPlane(-2 * hw)
vacLeft = openmc.XPlane(-2 * hw)
vacRight = openmc.XPlane(2 * hw)
vacFront = openmc.YPlane(-2 * hw)
vacBack = openmc.YPlane(2 * hw)

vacTop.boundary_type = 'vacuum'
vacBottom.boundary_type = 'vacuum'
vacLeft.boundary_type = 'vacuum'
vacRight.boundary_type = 'vacuum'
vacFront.boundary_type = 'vacuum'
vacBack.boundary_type = 'vacuum'

wpuCell = openmc.Cell()
wpuCell.region = -wpuSurface
wpuCell.fill = wpuMat

watCell = openmc.Cell()
watCell.region = +wpuSurface & +watBottom & -watTop & +watLeft & -watRight & +watFront & -watBack
watCell.fill = watMat

vacCell = openmc.Cell()
vacCell.region = +vacBottom & -vacTop & +vacLeft & -vacRight & +vacFront & -vacBack

root = openmc.Universe(cells = [wpuCell, watCell, vacCell])
geometry = openmc.Geometry(root)
geometry.export_to_xml(os.path.join(basepath, "geometry.xml"))

# Plot
if plot:
    print("Plotting geometry - can take a few seconds")

    plt.figure(figsize=(8, 8))
    root.plot(origin = (0, 0, 0),
              basis=('xz'),
              width=(2.3 * hw, 2.3* hw),
              pixels=(600, 600),
              seed = 1)
    plt.title("Water basin view (xz)")
    plt.savefig(os.path.join(basepath, "water-basin-xz.png"))
    plt.close()
    
###############################################################################
# Source & Settings
###############################################################################

sourcex = 0
sourcey = 0
sourcez = 0

bounds = [sourcex - wpuradius, sourcey - wpuradius, sourcez - wpuradius,
          sourcex + wpuradius, sourcey + wpuradius, sourcez + wpuradius]
uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable=True)
normalizationconstant = sum(pudf['wo'] * pudf['sfneutrons'])
source = []
for iso in pudf.index:
    tmpsrc = openmc.Source(space = uniform_dist, particle = 'neutron')
    tmpsrc.energy = openmc.stats.Watt(a=1/pudf.loc[iso, 'watt-a'], b=pudf.loc[iso, 'watt-b'])
    tmpsrc.strength = pudf.loc[iso, 'sfneutrons'] * pudf.loc[iso, 'wo'] / normalizationconstant
    source.append(tmpsrc)

sourcesum = 0
for iso in pudf.index:
    sourcesum += pudf.loc[iso, 'sfneutrons'] * pudf.loc[iso, 'wo']
print("A {:.1f} kg WPu weapon emits {:f} neutrons/s".format(wpumass / 1000, sourcesum * wpumass))
if time != 0:
    particles = int(sourcesum * wpumass * time)
    print("Simulation with {:n} particles/batch".format(particles))
    
settings = openmc.Settings()
settings.run_mode = 'fixed source'
settings.inactive = 0
settings.batches = batches
settings.particles = particles

if survival:
    settings.survival_biasing = True
settings.source = source

settings.export_to_xml(os.path.join(basepath, "settings.xml"))

###############################################################################
# Tallies
###############################################################################

tallies = openmc.Tallies()

xmin = -hw
xmax = +hw
ymin = -hw
ymax = +hw
zmin = hw
zmax = hw + 10

xstep = 10
tallyxmin = math.floor(xmin / xstep) * xstep
tallyxmax = math.ceil(xmax / xstep) * xstep 
ystep = 10
tallyymin = math.floor(ymin / ystep) * ystep
tallyymax = math.ceil(ymax / ystep) * ystep 
zstep = 10
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

meshfilter = openmc.MeshFilter(regularmesh)
tally = openmc.Tally(name = "flux, regular mesh")
tally.filters = [meshfilter]
tally.scores = ["flux"]
tallies.append(tally)

meshsurfacefilter = openmc.MeshSurfaceFilter(regularmesh)
tally = openmc.Tally(name = "current, regular mesh")
tally.filters = [meshsurfacefilter]
tally.scores = ["current"]
tallies.append(tally)

tallies.export_to_xml(os.path.join(basepath, "tallies.xml"))
