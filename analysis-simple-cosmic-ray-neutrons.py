import openmc
import json
import os
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib
import math
import copy
import sys

import helper

if len(sys.argv) == 2:
    basepath = sys.argv[1]
else:
    basepath = "cluster-results/sicorayn-20210726-10M"

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

puvec = helper.createwgpuvec()
age = weaponage
puvec.createagedvector([age])
sourcesum = 0
for iso in puvec.pudf.index:
    sourcesum += puvec.puagedf.loc[(iso, age), 'sfneutrons'] * puvec.puagedf.loc[(iso, age), 'wo']
nps = sourcesum * pumass
print("The WPu weapon source ({:.2f} years old, {:.1f}g) emits {:.1f} neutrons/s".format(age, pumass, nps))

################################################################################
# soil
soilsp = openmc.StatePoint(os.path.join(soilpath, "statepoint.{:d}.h5".format(batches)))
soilfettersp = openmc.StatePoint(os.path.join(soilfetterpath, "statepoint.{:d}.h5".format(batches)))

soiltally = soilsp.get_tally(name = 'flux, regular mesh')

cosmicsettings = openmc.Settings.from_xml(os.path.join(soilpath, "settings.xml"))
parameters = cosmicsettings.source[0].parameters
cosmicraywidth = float(parameters.split("subboxLength ")[-1].split(" ")[0])
cosmicneutrons = cosmicneutronsperm2 * cosmicraywidth ** 2

(xwidth, ywidth, zwidth) = soiltally.filters[0].mesh.dimension
soilfluxdata = soiltally.get_reshaped_data()[:,0,0]
soilfluxdata = soilfluxdata.reshape((zwidth, ywidth, xwidth))
soilfluxdata *= cosmicneutrons

plotcols = math.ceil(math.sqrt(len(soilfluxdata)))
plotrows = math.ceil(len(soilfluxdata) / plotcols)
fig, ax = plt.subplots(nrows = plotrows, ncols = plotcols, squeeze = False, figsize=(20, 20))

# thanks to https://stackoverflow.com/questions/9455044/problems-with-zeros-in-matplotlib-colors-lognorm
my_cmap = copy.copy(matplotlib.cm.get_cmap('viridis')) # copy the default cmap
my_cmap.set_bad(my_cmap(0))

vmin = soilfluxdata[soilfluxdata > 0].min()
vmax = soilfluxdata.max()

for i in range(plotcols * plotrows):
    rowidx = i // plotcols
    colidx = i % plotcols
    print(i, rowidx, colidx)
    if i < len(soilfluxdata):
        im = ax[rowidx, colidx].imshow(soilfluxdata[i], norm=LogNorm(vmin = vmin, vmax = vmax), cmap=my_cmap)
        ax[rowidx, colidx].title.set_text("z={:d}, Max: {:.3e}".format(i, soilfluxdata[i].max()))
    else:
        fig.delaxes(ax[rowidx, colidx])

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(im, cax=cbar_ax)
        
plt.savefig(os.path.join(plotpath, "soil-flux.png"))
plt.close()
#plt.show()

soiltally = soilsp.get_tally(name = 'current, regular mesh')
tempdf = soiltally.get_pandas_dataframe()
soilcurrentdf = tempdf[tempdf[('mesh 1', 'surf')].isin(['x-min in', 'x-max in', 'y-min in', 'y-max in', 'z-min in', 'z-max in'])]

soilcurrentdf = soilcurrentdf.groupby([('mesh 1', 'x'), ('mesh 1', 'y'), ('mesh 1', 'z')]).sum()
soilcurrentdf.reset_index(inplace = True)

soilcurrentdata = soilcurrentdf.sort_values([('mesh 1', 'z'), ('mesh 1', 'y'), ('mesh 1', 'x')])['mean'].to_numpy()
soilcurrentdata = soilcurrentdata.reshape((zwidth, ywidth, xwidth))
soilcurrentdata *= cosmicneutrons

vmin = 0.1
vmax = 10 ** math.ceil(math.log(cosmicneutronsperm2, 10))
print(vmin, vmax)

plotcols = math.ceil(math.sqrt(len(soilcurrentdata)))
plotrows = math.ceil(len(soilcurrentdata) / plotcols)
fig, ax = plt.subplots(nrows = plotrows, ncols = plotcols, squeeze = False, figsize=(20, 20))
for i in range(plotcols * plotrows):
    rowidx = i // plotcols
    colidx = i % plotcols
    print(i, rowidx, colidx)
    if i < len(soilcurrentdata):
        im = ax[rowidx, colidx].imshow(soilcurrentdata[i], norm=LogNorm(vmin = vmin, vmax = vmax), cmap=my_cmap)
        ax[rowidx, colidx].title.set_text("z={:d}, Max: {:.3e}".format(i, soilcurrentdata[i].max()))
    else:
        fig.delaxes(ax[rowidx, colidx])

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(im, cax=cbar_ax)
        

plt.savefig(os.path.join(plotpath, "soil-current.png"))
#plt.close()
plt.show()

soilfettertally = soilfettersp.get_tally(name = 'current, regular mesh')
tempdf = soilfettertally.get_pandas_dataframe()
soilfettercurrentdf = tempdf[tempdf[('mesh 1', 'surf')].isin(['x-min in', 'x-max in', 'y-min in', 'y-max in', 'z-min in', 'z-max in'])]

soilfettercurrentdf = soilfettercurrentdf.groupby([('mesh 1', 'x'), ('mesh 1', 'y'), ('mesh 1', 'z')]).sum()
soilfettercurrentdf.reset_index(inplace = True)

soilfettercurrentdata = soilfettercurrentdf.sort_values([('mesh 1', 'z'), ('mesh 1', 'y'), ('mesh 1', 'x')])['mean'].to_numpy()
soilfettercurrentdata = soilfettercurrentdata.reshape((zwidth, ywidth, xwidth))
soilfettercurrentdata *= nps

vmin = 0.1
vmax = 10 ** math.ceil(math.log(cosmicneutronsperm2, 10))
print(vmin, vmax)

plotcols = math.ceil(math.sqrt(len(soilfettercurrentdata)))
plotrows = math.ceil(len(soilfettercurrentdata) / plotcols)
fig, ax = plt.subplots(nrows = plotrows, ncols = plotcols, squeeze = False, figsize=(20, 20))
for i in range(plotcols * plotrows):
    rowidx = i // plotcols
    colidx = i % plotcols
    print(i, rowidx, colidx)
    if i < len(soilfettercurrentdata):
        im = ax[rowidx, colidx].imshow(soilfettercurrentdata[i], norm=LogNorm(vmin = vmin, vmax = vmax), cmap=my_cmap)
        ax[rowidx, colidx].title.set_text("z={:d}, Max: {:.3e}".format(i, soilfettercurrentdata[i].max()))
    else:
        fig.delaxes(ax[rowidx, colidx])
        
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(im, cax=cbar_ax)

plt.savefig(os.path.join(plotpath, "soil-and-fetter-current.png"))
#plt.close()
plt.show()

z = 6
fig, ax = plt.subplots(nrows = 1, ncols = 3, squeeze = False, figsize=(21, 7))
fig.subplots_adjust(left=0.02, bottom=0.06, right=0.95, top=0.94, wspace=0.05)

im = ax[0, 0].imshow(soilcurrentdata[z], norm=LogNorm(vmin = vmin, vmax = vmax), cmap=my_cmap)
ax[0, 0].title.set_text("Cosmic ray neutron background level")
fig.colorbar(im, ax = ax[0, 0])

im = ax[0, 1].imshow(soilfettercurrentdata[z], norm=LogNorm(vmin = vmin, vmax = vmax), cmap=my_cmap)
ax[0, 1].title.set_text("Neutron level from nuclear weapon")
fig.colorbar(im, ax = ax[0, 1])

tdata = soilcurrentdata[z] / (soilfettercurrentdata[z] ** 2)
im = ax[0, 2].imshow(tdata, norm=LogNorm(), cmap=my_cmap)
ax[0, 2].title.set_text("One sigma measurement time (s)")
fig.colorbar(im, ax = ax[0, 2])

plt.savefig(os.path.join(plotpath, "soil-measurement-times.png"))
plt.show()
