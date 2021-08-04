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
    basepath = "cluster-results/jovanovic-erickson-5000-1seconds"

rate = 311324 # neutrons per second in reality

with open(os.path.join(basepath, "calculation.json"), 'r') as f:
    calcsettings = json.load(f)
    batches = calcsettings['batches']
    f.close()

sp = openmc.StatePoint(os.path.join(basepath, "statepoint.{:d}.h5".format(batches)))
tally = sp.get_tally(name = 'flux, regular mesh')
#flux = tally.get_slice(scores=['flux'])

(xwidth, ywidth, zwidth) = tally.filters[0].mesh.dimension
fluxdata = tally.get_reshaped_data()[:,0,0]
fluxdata = fluxdata.reshape((zwidth, ywidth, xwidth))
fluxdata = fluxdata * rate

fig, ax = plt.subplots(1, 1, figsize=(20, 20))

# thanks to https://stackoverflow.com/questions/9455044/problems-with-zeros-in-matplotlib-colors-lognorm
my_cmap = copy.copy(matplotlib.cm.get_cmap('viridis')) # copy the default cmap
my_cmap.set_bad("black")

i = 0
#im = ax.imshow(fluxdata[i], norm=LogNorm(vmin=1e-10, vmax=1e-4), cmap=my_cmap)
im = ax.imshow(fluxdata[i], norm=LogNorm(), cmap=my_cmap)
ax.title.set_text("z={:d}, Max: {:.3e}".format(i, fluxdata[i].max()))
fig.colorbar(im, ax = ax)
plt.savefig("weapon-in-water-neutron-flux.png")
#plt.close()
plt.show()

tally = sp.get_tally(name = 'current, regular mesh')
tempdf = tally.get_pandas_dataframe()
currentdf = tempdf[tempdf[('mesh 1', 'surf')] == 'z-min in']
#currentoutdf = tempdf[tempdf[('mesh 2', 'surf')].isin(['x-min out', 'x-max out', 'y-min out', 'y-max out', 'z-min out', 'z-max out'])]

# currentdf = currentdf.groupby([('mesh 2', 'x'), ('mesh 2', 'y'), ('mesh 2', 'z')]).sum()
# currentdf.reset_index(inplace = True)

currentdata = currentdf.sort_values([('mesh 1', 'z'), ('mesh 1', 'y'), ('mesh 1', 'x')])['mean'].to_numpy()
currentdata = currentdata.reshape((zwidth, ywidth, xwidth))
currentdata = currentdata * rate
currentdata = currentdata / 100 # 10 by 10 cm surfaces

# i = 4
# fig, ax = plt.subplots(figsize=(5, 5))
# ax.imshow(fluxdata[i], norm=LogNorm(vmin=0.000001, vmax=10), cmap=my_cmap)
# ax.title.set_text("z={:d}, Max: {:.3e}".format(i, fluxdata[i].max()))
# plt.savefig("buechel-storage-neutron-flux-4m.png")    
# plt.show()

fig, ax = plt.subplots(1, 1, figsize=(20, 20))

im = ax.imshow(currentdata[i], norm=LogNorm(), cmap=my_cmap)
ax.title.set_text("z={:d}, Max: {:.3e}".format(i, currentdata[i].max()))

cbar = fig.colorbar(im, ax = ax)
cbar.set_label("neutrons / cm^2 / s")

plt.savefig("weapon-in-water-neutron-current.png")
#plt.close()
plt.show()

#df2 = df[df[('mesh 2', 'surf')].isin(['x-min in', 'x-max in', 'y-min in', 'y-max in', 'z-min in', 'z-max in'])]
