import openmc
import json
import os
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib
import math
import copy
import sys
import numpy as np
import pandas as pd

import helper

if len(sys.argv) == 2:
    basepath = sys.argv[1]
else:
    basepath = "cluster-results/sicorayn-20210727-20M/"

concretepath = os.path.join(basepath, "concrete")
waterpath = os.path.join(basepath, "water")
soilpath = os.path.join(basepath, "soil")
concrete2100path = os.path.join(basepath, "concrete-2100")
water2100path = os.path.join(basepath, "water-2100")
soil2100path = os.path.join(basepath, "soil-2100")
concretefetterpath = os.path.join(basepath, "concrete-fetter")
waterfetterpath = os.path.join(basepath, "water-fetter")
soilfetterpath = os.path.join(basepath, "soil-fetter")

plotpath = os.path.join("plots", "sicorayn")
if not os.path.exists(plotpath):
    os.mkdir(plotpath)

batches = 10
cosmicneutronsperm2 = 120
pumass = 4000
weaponage = 0
weaponz = 6

puvec = helper.createwgpuvec()
age = weaponage
puvec.createagedvector([age])
sourcesum = 0
for iso in puvec.pudf.index:
    sourcesum += puvec.puagedf.loc[(iso, age), 'sfneutrons'] * puvec.puagedf.loc[(iso, age), 'wo']
nps = sourcesum * pumass
print("The WPu weapon source ({:.2f} years old, {:.1f}g) emits {:.1f} neutrons/s".format(age, pumass, nps))

n_surfs = 12
stride = 1
def _repeat_and_tile(bins, repeat_factor, data_size):
    filter_bins = np.repeat(bins, repeat_factor)
    tile_factor = data_size // len(filter_bins)
    return np.tile(filter_bins, tile_factor)

def meshsurfacetodataframe(tally, xwidth, ywidth, zwidth, factor = 1):
    data = tally.get_reshaped_data().flatten()
    data_size = len(data) // n_surfs
    data = data.reshape((data_size, n_surfs)).transpose()
    filter_dict = {}
    filter_dict['x'] = _repeat_and_tile(np.arange(xwidth), stride, data_size)
    filter_dict['y'] = _repeat_and_tile(np.arange(ywidth), xwidth * stride, data_size)
    filter_dict['z'] = _repeat_and_tile(np.arange(zwidth), xwidth * ywidth * stride, data_size)
    filterdf = pd.DataFrame(filter_dict)

    for i in range(n_surfs):
        name = openmc.filter._CURRENT_NAMES[i]
        filterdf[name] = data[i] * factor

    filterdf['in'] = filterdf[['x-min in', 'x-max in', 'y-min in', 'y-max in', 'z-min in', 'z-max in']].sum(axis = 1)
    filterdf['out'] = filterdf[['x-min out', 'x-max out', 'y-min out', 'y-max out', 'z-min out', 'z-max out']].sum(axis = 1)
    filterdf['xyin'] = filterdf[['x-min in', 'x-max in', 'y-min in', 'y-max in']].sum(axis = 1)
    filterdf['xyout'] = filterdf[['x-min out', 'x-max out', 'y-min out', 'y-max out']].sum(axis = 1)

    return filterdf

# thanks to https://stackoverflow.com/questions/9455044/problems-with-zeros-in-matplotlib-colors-lognorm
my_cmap = copy.copy(matplotlib.cm.get_cmap('viridis')) # copy the default cmap
my_cmap.set_bad(my_cmap(0))

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

centerlines = {
    "source": [],
    "z": [],
    "data": []
}


################################################################################
print("LOADING SOIL DATA")

# load statepoint
sp = openmc.StatePoint(os.path.join(soilpath, "statepoint.{:d}.h5".format(batches)))
fsp = openmc.StatePoint(os.path.join(soilfetterpath, "statepoint.{:d}.h5".format(batches)))

settings = openmc.Settings.from_xml(os.path.join(soilpath, "settings.xml"))
parameters = settings.source[0].parameters
soilwidth = float(parameters.split("subboxLength ")[-1].split(" ")[0])
soilcosmicneutrons = cosmicneutronsperm2 * soilwidth ** 2

# load meshsurface tally for cosmic rays (current)
print("mesh surface tally (current)")
tally = sp.get_tally(name = 'current, regular mesh')
(xwidth, ywidth, zwidth) = tally.filters[0].mesh.dimension
soilcosmicdf = meshsurfacetodataframe(tally, xwidth, ywidth, zwidth, soilcosmicneutrons)
soilcosmicdata = soilcosmicdf['in'].values.reshape(zwidth, ywidth, xwidth)

ycenter = math.floor(ywidth / 2)
xcenter = math.floor(xwidth / 2)
for z in range(zwidth):
    centerlines['z'].append(z)
    centerlines['source'].append("soil")
    centerlines['data'].append(soilcosmicdf.query("z=={} & y=={}".format(z, ycenter))['in'].values)
    centerlines['z'].append(z)
    centerlines['source'].append("soil-shielded")
    centerlines['data'].append(soilcosmicdf.query("z=={} & y=={}".format(z, ycenter))['xyin'].values)

# # use this to compare / debug
# comparedf = meshsurfacetodataframe(tally, xwidth, ywidth, zwidth)
# tempdf = tally.get_pandas_dataframe()
# t1 = tempdf[tempdf[('mesh 1', 'surf')] == 'x-min out']['mean'].values
# t2 = comparedf['x-min out'].values
# print(t1)
# print(t2)

tally = fsp.get_tally(name = 'current, regular mesh')
(xwidth, ywidth, zwidth) = tally.filters[0].mesh.dimension
soilfetterdf = meshsurfacetodataframe(tally, xwidth, ywidth, zwidth, nps)
soilfetterdata = soilfetterdf['in'].values.reshape(zwidth, ywidth, xwidth)

ycenter = math.floor(ywidth / 2)
xcenter = math.floor(xwidth / 2)
for z in range(zwidth):
    centerlines['z'].append(z)
    centerlines['source'].append("fetter-on-soil")
    centerlines['data'].append(soilfetterdf.query("z=={} & y=={}".format(z, ycenter))['in'].values)
    centerlines['z'].append(z)
    centerlines['source'].append("fetter-on-soil-shielded")
    centerlines['data'].append(soilfetterdf.query("z=={} & y=={}".format(z, ycenter))['xyin'].values)

weaponcellout = soilfetterdf.query("z=={} & y=={} & x=={}".format(weaponz, ycenter, xcenter))['out'].values[0]
weaponcellin = soilfetterdf.query("z=={} & y=={} & x=={}".format(weaponz, ycenter, xcenter))['in'].values[0]
soilmax = weaponcellout - weaponcellin

print("The weapon was simulated in mesh cell x={}, y={}, z={} (including 0 in count).".format(xcenter, ycenter, weaponz))
print("Outgoing current {:e}".format(weaponcellout))
print("Incoming current {:e}".format(weaponcellin))
print("Net current {:e}".format(weaponcellout - weaponcellin))

centerlinedf = pd.DataFrame(centerlines)
centerlinedf.set_index(['source', 'z'], inplace=True)

################################################################################
print("LOADING CONCRETE DATA")

# load statepoint
sp = openmc.StatePoint(os.path.join(concretepath, "statepoint.{:d}.h5".format(batches)))
fsp = openmc.StatePoint(os.path.join(concretefetterpath, "statepoint.{:d}.h5".format(batches)))

settings = openmc.Settings.from_xml(os.path.join(concretepath, "settings.xml"))
parameters = settings.source[0].parameters
concretewidth = float(parameters.split("subboxLength ")[-1].split(" ")[0])
concretecosmicneutrons = cosmicneutronsperm2 * concretewidth ** 2

# load meshsurface tally for cosmic rays (current)
print("mesh surface tally (current)")
tally = sp.get_tally(name = 'current, regular mesh')
(xwidth, ywidth, zwidth) = tally.filters[0].mesh.dimension
concretecosmicdf = meshsurfacetodataframe(tally, xwidth, ywidth, zwidth, concretecosmicneutrons)
concretecosmicdata = concretecosmicdf['in'].values.reshape(zwidth, ywidth, xwidth)

ycenter = math.floor(ywidth / 2)
xcenter = math.floor(xwidth / 2)
for z in range(zwidth):
    centerlines['z'].append(z)
    centerlines['source'].append("concrete")
    centerlines['data'].append(concretecosmicdf.query("z=={} & y=={}".format(z, ycenter))['in'].values)
    centerlines['z'].append(z)
    centerlines['source'].append("concrete-shielded")
    centerlines['data'].append(concretecosmicdf.query("z=={} & y=={}".format(z, ycenter))['xyin'].values)

# # use this to compare / debug
# comparedf = meshsurfacetodataframe(tally, xwidth, ywidth, zwidth)
# tempdf = tally.get_pandas_dataframe()
# t1 = tempdf[tempdf[('mesh 1', 'surf')] == 'x-min in']['mean'].values
# t2 = comparedf['x-min in'].values
# print(t1)
# print(t2)

tally = fsp.get_tally(name = 'current, regular mesh')
(xwidth, ywidth, zwidth) = tally.filters[0].mesh.dimension
concretefetterdf = meshsurfacetodataframe(tally, xwidth, ywidth, zwidth, nps)
concretefetterdata = concretefetterdf['in'].values.reshape(zwidth, ywidth, xwidth)

ycenter = math.floor(ywidth / 2)
xcenter = math.floor(xwidth / 2)
for z in range(zwidth):
    centerlines['z'].append(z)
    centerlines['source'].append("fetter-on-concrete")
    centerlines['data'].append(concretefetterdf.query("z=={} & y=={}".format(z, ycenter))['in'].values)
    centerlines['z'].append(z)
    centerlines['source'].append("fetter-on-concrete-shielded")
    centerlines['data'].append(concretefetterdf.query("z=={} & y=={}".format(z, ycenter))['xyin'].values)

weaponcellout = concretefetterdf.query("z=={} & y=={} & x=={}".format(weaponz, ycenter, xcenter))['out'].values[0]
weaponcellin = concretefetterdf.query("z=={} & y=={} & x=={}".format(weaponz, ycenter, xcenter))['in'].values[0]
concretemax = weaponcellout - weaponcellin

print("The weapon was simulated in mesh cell x={}, y={}, z={} (including 0 in count).".format(xcenter, ycenter, weaponz))
print("Outgoing current {:e}".format(weaponcellout))
print("Incoming current {:e}".format(weaponcellin))
print("Net current {:e}".format(weaponcellout - weaponcellin))

centerlinedf = pd.DataFrame(centerlines)
centerlinedf.set_index(['source', 'z'], inplace=True)

################################################################################
print("LOADING WATER DATA")

# load statepoint
sp = openmc.StatePoint(os.path.join(waterpath, "statepoint.{:d}.h5".format(batches)))
fsp = openmc.StatePoint(os.path.join(waterfetterpath, "statepoint.{:d}.h5".format(batches)))

settings = openmc.Settings.from_xml(os.path.join(waterpath, "settings.xml"))
parameters = settings.source[0].parameters
waterwidth = float(parameters.split("subboxLength ")[-1].split(" ")[0])
watercosmicneutrons = cosmicneutronsperm2 * waterwidth ** 2

# load meshsurface tally for cosmic rays (current)
print("mesh surface tally (current)")
tally = sp.get_tally(name = 'current, regular mesh')
(xwidth, ywidth, zwidth) = tally.filters[0].mesh.dimension
watercosmicdf = meshsurfacetodataframe(tally, xwidth, ywidth, zwidth, watercosmicneutrons)
watercosmicdata = watercosmicdf['in'].values.reshape(zwidth, ywidth, xwidth)

ycenter = math.floor(ywidth / 2)
xcenter = math.floor(xwidth / 2)
for z in range(zwidth):
    centerlines['z'].append(z)
    centerlines['source'].append("water")
    centerlines['data'].append(watercosmicdf.query("z=={} & y=={}".format(z, ycenter))['in'].values)
    centerlines['z'].append(z)
    centerlines['source'].append("water-shielded")
    centerlines['data'].append(watercosmicdf.query("z=={} & y=={}".format(z, ycenter))['xyin'].values)

# # use this to compare / debug
# comparedf = meshsurfacetodataframe(tally, xwidth, ywidth, zwidth)
# tempdf = tally.get_pandas_dataframe()
# t1 = tempdf[tempdf[('mesh 1', 'surf')] == 'x-min out']['mean'].values
# t2 = comparedf['x-min out'].values
# print(t1)
# print(t2)

tally = fsp.get_tally(name = 'current, regular mesh')
(xwidth, ywidth, zwidth) = tally.filters[0].mesh.dimension
waterfetterdf = meshsurfacetodataframe(tally, xwidth, ywidth, zwidth, nps)
waterfetterdata = waterfetterdf['in'].values.reshape(zwidth, ywidth, xwidth)

ycenter = math.floor(ywidth / 2)
xcenter = math.floor(xwidth / 2)
for z in range(zwidth):
    centerlines['z'].append(z)
    centerlines['source'].append("fetter-on-water")
    centerlines['data'].append(waterfetterdf.query("z=={} & y=={}".format(z, ycenter))['in'].values)
    centerlines['z'].append(z)
    centerlines['source'].append("fetter-on-water-shielded")
    centerlines['data'].append(waterfetterdf.query("z=={} & y=={}".format(z, ycenter))['xyin'].values)

weaponcellout = waterfetterdf.query("z=={} & y=={} & x=={}".format(weaponz, ycenter, xcenter))['out'].values[0]
weaponcellin = waterfetterdf.query("z=={} & y=={} & x=={}".format(weaponz, ycenter, xcenter))['in'].values[0]
watermax = weaponcellout - weaponcellin

print("The weapon was simulated in mesh cell x={}, y={}, z={} (including 0 in count).".format(xcenter, ycenter, weaponz))
print("Outgoing current {:e}".format(weaponcellout))
print("Incoming current {:e}".format(weaponcellin))
print("Net current {:e}".format(weaponcellout - weaponcellin))

centerlinedf = pd.DataFrame(centerlines)
centerlinedf.set_index(['source', 'z'], inplace=True)

#*******************************************************************************
# make a comprehensive plot
# left: weapon / background
# center: shielded / unshielded
# right: measurement time

factor = 3
fig, ax = plt.subplots(nrows = 1, ncols = 3, squeeze = False, figsize=(3 * factor, 1 * factor))

xvals = np.arange(xwidth) - xcenter
ax[0][0].plot(xvals, centerlinedf.loc[('soil', weaponz), 'data'], label="Soil", color=colors[0])
ax[0][0].plot(xvals, centerlinedf.loc[('fetter-on-soil', weaponz), 'data'], "--", color=colors[0])

ax[0][0].plot(xvals, centerlinedf.loc[('concrete', weaponz), 'data'], label="Concrete", color=colors[1])
ax[0][0].plot(xvals, centerlinedf.loc[('fetter-on-concrete', weaponz), 'data'], "--", color=colors[1])

ax[0][0].plot(xvals, centerlinedf.loc[('water', weaponz), 'data'], label="Water", color=colors[2])
ax[0][0].plot(xvals, centerlinedf.loc[('fetter-on-water', weaponz), 'data'], "--", color=colors[2])

ax[0][0].set_ylabel("Current into 1 m$^3$ [neutrons/s]")
ax[0][0].set_xlabel("Horizontal Distance [m]")

ax[0][0].set_ylim(-100, 1100)
ax[0][0].set_xlim(0, 100)
ax[0][0].grid()
ax[0][0].legend()

ax[0][0].title.set_text("Different Surface Materials")

#*******************************************************************************
# plot difference shielded / unshielded
ax[0][1].plot(xvals, centerlinedf.loc[('soil', weaponz), 'data'], label="Unshielded", color=colors[0])
ax[0][1].plot(xvals, centerlinedf.loc[('fetter-on-soil', weaponz), 'data'], "--", color=colors[0])

ax[0][1].plot(xvals, centerlinedf.loc[('soil-shielded', weaponz), 'data'], label="Shielded", color=colors[1])
ax[0][1].plot(xvals, centerlinedf.loc[('fetter-on-soil-shielded', weaponz), 'data'], "--", color=colors[1])

ax[0][1].set_ylabel("Current into 1 m$^3$ [neutrons/s]")
ax[0][1].set_xlabel("Horizontal Distance [m]")

ax[0][1].set_ylim(-100, 1100)
ax[0][1].set_xlim(0, 100)
ax[0][1].grid()
ax[0][1].legend()

# ax[0][1].title.set_text("Detector with/without z-direction Shielding")
ax[0][1].title.set_text("Effect of z-direction Shielding")

#*******************************************************************************
# plot measurement time - five sigma

m = 5
tdata1 = m ** 2 * centerlinedf.loc[('soil', weaponz), 'data'] / (centerlinedf.loc[('fetter-on-soil', weaponz), 'data'] ** 2)
tdata1[xcenter] = np.nan # no useful calculation possible here
ax[0][2].plot(xvals, tdata1, label = "Soil - Unshielded")

tdata2 = m ** 2 * centerlinedf.loc[('soil-shielded', weaponz), 'data'] / (centerlinedf.loc[('fetter-on-soil-shielded', weaponz), 'data'] ** 2)
tdata2[xcenter] = np.nan # no useful calculation possible here
ax[0][2].plot(xvals, tdata2, label = "Soil - Shielded")

ax[0][2].plot(xvals, tdata1 - tdata2, label = "Difference")

ax[0][2].set_ylabel("Time [s]")
ax[0][2].set_xlabel("Horizontal Distance [m]")

ax[0][2].axhline(60, color = colors[6], linewidth = 2)
ax[0][2].axhline(3600, color = colors[6], linewidth = 2)
ax[0][2].axhline(3600 * 24, color = colors[6], linewidth = 2)

ax[0][2].set_xlim(0, 100)
ax[0][2].set_yscale("log")
#ax[0][2].set_ylim(0, 100)
ax[0][2].grid()
ax[0][2].legend()

ax[0][2].title.set_text("Measurement Time")

fig.tight_layout()

plt.savefig(os.path.join("plots-for-paper", "simple-cosmic-centerline.pdf"))
plt.savefig(os.path.join("plots-for-paper", "simple-cosmic-centerline.png"))

plt.show()

#*******************************************************************************
# plot cosmic / fetter on soil for presentation
matplotlib.rcParams['font.sans-serif'] = "Fira Sans"
matplotlib.rcParams['font.family'] = "sans-serif"

factor = 2.7

fig, ax = plt.subplots(nrows = 1, ncols = 2, squeeze = False, figsize=(2.4 * factor, 1 * factor))
fig.subplots_adjust(left=0.02, bottom=0.06, right=0.95, top=0.94, wspace=0.05)

vmin = soilcosmicdata[soilcosmicdata > 0].min()
vmax = 10 ** math.ceil(math.log(soilcosmicdata.max(), 10) + 2)
vmin = 1e-1
vmax = 1e4
print(vmax, vmin)

im = ax[0, 0].imshow(soilcosmicdata[weaponz], norm=LogNorm(vmin = vmin, vmax = vmax), cmap=my_cmap)
ax[0, 0].title.set_text("Cosmic-ray neutron background")
ax[0, 0].set_ylabel("[meters]")
ax[0, 0].set_xlabel("[meters]")
fig.colorbar(im, ax = ax[0, 0])

im = ax[0, 1].imshow(soilfetterdata[weaponz], norm=LogNorm(vmin = vmin, vmax = vmax), cmap=my_cmap)
ax[0, 1].title.set_text("Signal from nuclear weapon")
ax[0, 1].set_xlabel("[meters]")
fig.colorbar(im, ax = ax[0, 1], label = "Current into 1 m$^3$ [neutrons/s]")#"neutrons entering 1m$^3$/s")

fig.tight_layout()
plt.savefig(os.path.join("plots-for-presentation", "soil.pdf"), transparent = True)
plt.savefig(os.path.join("plots-for-presentation", "soil.png"), transparent = True)
plt.show()

#*******************************************************************************
# make a comprehensive plot for presentation
# left: weapon / background
# center: shielded / unshielded
# right: measurement time

matplotlib.rcParams['font.sans-serif'] = "Fira Sans"
matplotlib.rcParams['font.family'] = "sans-serif"

factor = 3.2
fig, ax = plt.subplots(nrows = 1, ncols = 2, squeeze = False, figsize=(2 * factor, 1 * factor))

xvals = np.arange(xwidth) - xcenter
ax[0][0].plot(xvals, centerlinedf.loc[('soil', weaponz), 'data'], label="Soil", color=colors[0])
ax[0][0].plot(xvals, centerlinedf.loc[('soil-shielded', weaponz), 'data'], "--", label="Soil - Shielded", color=colors[0])

ax[0][0].plot(xvals, centerlinedf.loc[('concrete', weaponz), 'data'], label="Concrete", color=colors[1])

ax[0][0].plot(xvals, centerlinedf.loc[('water', weaponz), 'data'], label="Water", color=colors[2])

ax[0][0].plot(xvals, centerlinedf.loc[('fetter-on-soil', weaponz), 'data'], label = "Weapon Signal", color=colors[4])

ax[0][0].set_ylabel("Current into 1 m$^3$ [neutrons/s]")
ax[0][0].set_xlabel("Horizontal Distance [m]")

ax[0][0].set_ylim(-100, 1100)
ax[0][0].set_xlim(0, 100)
ax[0][0].grid()
ax[0][0].legend(ncol=1)

ax[0][0].title.set_text("Different Surface Materials")

#*******************************************************************************
# plot measurement time - five sigma

m = 5
tdata1 = m ** 2 * centerlinedf.loc[('soil', weaponz), 'data'] / (centerlinedf.loc[('fetter-on-soil', weaponz), 'data'] ** 2)
tdata1[xcenter] = np.nan # no useful calculation possible here
ax[0][1].plot(xvals, tdata1, label = "Soil - Unshielded")

tdata2 = m ** 2 * centerlinedf.loc[('soil-shielded', weaponz), 'data'] / (centerlinedf.loc[('fetter-on-soil-shielded', weaponz), 'data'] ** 2)
tdata2[xcenter] = np.nan # no useful calculation possible here
ax[0][1].plot(xvals, tdata2, label = "Soil - Shielded")

ax[0][1].plot(xvals, tdata1 - tdata2, label = "Difference")

ax[0][1].set_ylabel("Time [s]")
ax[0][1].set_xlabel("Horizontal Distance [m]")

ax[0][1].axhline(60, color = colors[6], linewidth = 2)
ax[0][1].axhline(3600, color = colors[6], linewidth = 2)
ax[0][1].axhline(3600 * 24, color = colors[6], linewidth = 2)

ax[0][1].set_xlim(0, 100)
ax[0][1].set_yscale("log")
#ax[0][1].set_ylim(0, 100)
ax[0][1].grid()
ax[0][1].legend()

ax[0][1].title.set_text("Measurement Time")

fig.tight_layout()

plt.savefig(os.path.join("plots-for-presentation", "simple-cosmic-centerline.pdf"), transparent = True)
plt.savefig(os.path.join("plots-for-presentation", "simple-cosmic-centerline.png"), transparent = True)

plt.show()


