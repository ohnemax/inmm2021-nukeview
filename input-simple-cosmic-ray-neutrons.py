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
basepath = "simple-cosmic-ray-neutrons"
particles = 10000
batches = 10
plot = False
survival = False
discard = False
maxenergy = 2e7
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
parser.add_argument("-m", "--maxenergy", type = float,
                    help="source particles cutoff energy?")
parser.add_argument("-d", "--discard", action="store_true",
                    help="discard source particles above cutoff energy?")
parser.add_argument("-x", "--xsection",
                    help="specify cross section path")

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
if not args.maxenergy is None:
    maxenergy = args.maxenergy
if args.discard:
    discard = True
if not args.xsection is None:
    cspath = args.xsection

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

materiallist = [airMat, soilMat]

materials = openmc.Materials(materiallist)

# Check for cross section availability
lib = openmc.data.DataLibrary.from_xml(cspath)
for mat in materials:
    for nuclide in mat.nuclides:
        if lib.get_by_material(nuclide.name) is None:
            print("WARNING: Could not find {:s} in cross section library".format(nuclide.name))
            print("  The material contains {:e} {:s} of that nuclide.".format(nuclide.percent, nuclide.percent_type))
            print("  Will remove the nuclide from the material")
            mat.remove_nuclide(nuclide.name)
materials.cross_sections = cspath
materials.export_to_xml(os.path.join(basepath, "materials.xml"))

###############################################################################
# Geometry
###############################################################################

hw = 500
soilTop = openmc.ZPlane(0)
soilBottom = openmc.ZPlane(-hw)
soilLeft = openmc.XPlane(-hw)
soilRight = openmc.XPlane(hw)
soilFront = openmc.YPlane(-hw)
soilBack = openmc.YPlane(hw)

airTop = openmc.ZPlane(hw)
airBottom = openmc.ZPlane(0)
airLeft = openmc.XPlane(-hw)
airRight = openmc.XPlane(hw)
airFront = openmc.YPlane(-hw)
airBack = openmc.YPlane(hw)

soilBottom.boundary_type = 'vacuum'

soilLeft.boundary_type = 'vacuum'
soilRight.boundary_type = 'vacuum'
soilFront.boundary_type = 'vacuum'
soilBack.boundary_type = 'vacuum'

airLeft.boundary_type = 'vacuum'
airRight.boundary_type = 'vacuum'
airFront.boundary_type = 'vacuum'
airBack.boundary_type = 'vacuum'

airTop.boundary_type = 'vacuum'

soilCell = openmc.Cell()
soilCell.region = +soilBottom & -soilTop & +soilLeft & -soilRight & +soilFront & -soilBack
soilCell.fill = soilMat

airCell = openmc.Cell()
airCell.region = +airBottom & -airTop & +airLeft & -airRight & +airFront & -airBack
airCell.fill = airMat

root = openmc.Universe(cells = [soilCell, airCell])
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
    plt.title("Soil and air (xz)")
    plt.savefig(os.path.join(basepath, "soil-air-xz.png"))
    plt.close()
    
###############################################################################
# Source & Settings
###############################################################################
scriptdir = os.path.dirname(os.path.realpath(__file__))
if os.path.islink(os.path.join(basepath, "libsource.so")):
    os.unlink(os.path.join(basepath, "libsource.so"))
os.symlink(os.path.join(scriptdir, "../cry-with-openmc/build/libsource.so"), os.path.join(basepath, "libsource.so"))
if os.path.islink(os.path.join(basepath, "data")):
    os.unlink(os.path.join(basepath, "data"))
os.symlink(os.path.join(scriptdir, "../cry-with-openmc/build/CRY-1.7-prefix/src/CRY-1.7/data"), os.path.join(basepath, "data"))
source = openmc.Source(library = "./libsource.so")
if discard:
    discardstring = "discard"
else:
    discardstring = "limit"
latitude = 50.173833 # buechel airbase
source.parameters = "{:f} {:s} {:f} {:f} {:f} returnNeutrons 1 returnProtons 0 returnGammas 0 returnMuons 0 returnElectrons 0 returnPions 0 date 1-1-2008 latitude {:f} altitude 0 subboxLength {:f}".format(maxenergy, discardstring, 0, 0, hw - 1, latitude, 2 * hw / 100)
    
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
zmin = 0
zmax = hw

xstep = 100
tallyxmin = math.floor(xmin / xstep) * xstep
tallyxmax = math.ceil(xmax / xstep) * xstep 
ystep = 100
tallyymin = math.floor(ymin / ystep) * ystep
tallyymax = math.ceil(ymax / ystep) * ystep 
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

energyfilter = openmc.EnergyFilter(np.logspace(-3, 8, num=100))
meshfilter = openmc.MeshFilter(regularmesh)
tally = openmc.Tally(name = "flux, regular mesh")
tally.filters = [meshfilter, energyfilter]
tally.scores = ["flux"]
tallies.append(tally)

meshsurfacefilter = openmc.MeshSurfaceFilter(regularmesh)
tally = openmc.Tally(name = "current, regular mesh")
tally.filters = [meshsurfacefilter]
tally.scores = ["current"]
tallies.append(tally)

tallies.export_to_xml(os.path.join(basepath, "tallies.xml"))
