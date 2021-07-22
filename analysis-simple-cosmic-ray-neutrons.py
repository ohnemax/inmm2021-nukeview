import openmc
import json
import os
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib
import math
import copy
import sys

if len(sys.argv) == 2:
    basepath = sys.argv[1]
else:
    basepath = "simple-cosmic-ray-neutrons"

with open(os.path.join(basepath, "calculation.json"), 'r') as f:
    calcsettings = json.load(f)
    batches = calcsettings['batches']
    f.close()

sp = openmc.StatePoint(os.path.join(basepath, "statepoint.{:d}.h5".format(batches)))
tally = sp.get_tally(name = 'flux, regular mesh')
tally = tally.summation(filter_type = openmc.EnergyFilter, remove_filter=True)
#flux = tally.get_slice(scores=['flux'])

(xwidth, ywidth, zwidth) = tally.filters[0].mesh.dimension
fluxdata = tally.get_reshaped_data()[:,0,0]
fluxdata = fluxdata.reshape((zwidth, ywidth, xwidth))

plotcols = math.ceil(math.sqrt(len(fluxdata)))
plotrows = math.ceil(len(fluxdata) / plotcols)
fig, ax = plt.subplots(nrows = plotrows, ncols = plotcols, squeeze = False, figsize=(20, 20))

# thanks to https://stackoverflow.com/questions/9455044/problems-with-zeros-in-matplotlib-colors-lognorm
my_cmap = copy.copy(matplotlib.cm.get_cmap('viridis')) # copy the default cmap
my_cmap.set_bad(my_cmap(0))

for i in range(plotcols * plotrows):
    rowidx = i // plotcols
    colidx = i % plotcols
    print(i, rowidx, colidx)
    if i < len(fluxdata):
        ax[rowidx, colidx].imshow(fluxdata[i], norm=LogNorm(), cmap=my_cmap)
        ax[rowidx, colidx].title.set_text("z={:d}, Max: {:.3e}".format(i, fluxdata[i].max()))
    else:
        fig.delaxes(ax[rowidx, colidx])

plt.savefig("simple-cosmic-ray-neutrons-flux.png")
#plt.close()
plt.show()

tally = sp.get_tally(name = 'current, regular mesh')
tempdf = tally.get_pandas_dataframe()
currentdf = tempdf[tempdf[('mesh 1', 'surf')].isin(['x-min in', 'x-max in', 'y-min in', 'y-max in', 'z-min in', 'z-max in'])]

currentdf = currentdf.groupby([('mesh 1', 'x'), ('mesh 1', 'y'), ('mesh 1', 'z')]).sum()
currentdf.reset_index(inplace = True)

currentdata = currentdf.sort_values([('mesh 1', 'z'), ('mesh 1', 'y'), ('mesh 1', 'x')])['mean'].to_numpy()
currentdata = currentdata.reshape((zwidth, ywidth, xwidth))

plotcols = math.ceil(math.sqrt(len(currentdata)))
plotrows = math.ceil(len(currentdata) / plotcols)
fig, ax = plt.subplots(nrows = plotrows, ncols = plotcols, squeeze = False, figsize=(20, 20))
for i in range(plotcols * plotrows):
    rowidx = i // plotcols
    colidx = i % plotcols
    print(i, rowidx, colidx)
    if i < len(currentdata):
        ax[rowidx, colidx].imshow(currentdata[i], norm=LogNorm(), cmap=my_cmap)
        ax[rowidx, colidx].title.set_text("z={:d}, Max: {:.3e}".format(i, currentdata[i].max()))
    else:
        fig.delaxes(ax[rowidx, colidx])
        

plt.savefig("simple-cosmic-ray-neutrons-current.png")
#plt.close()
plt.show()

tally = sp.get_tally(name = 'flux, regular mesh')
tally = tally.summation(filter_type = openmc.MeshFilter, remove_filter=True)

fluxdf = tally.get_pandas_dataframe()
fluxdf.plot(x = 'energy low [eV]', y = 'mean',
            logx = True,
            logy = True)
plt.show()

#df2 = df[df[('mesh 2', 'surf')].isin(['x-min in', 'x-max in', 'y-min in', 'y-max in', 'z-min in', 'z-max in'])]
