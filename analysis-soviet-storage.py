import openmc
import json
import os
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib
import math
import copy

basepath = "sost-20210626-100k"
basepath = "cluster-results/sost-20210626-run-1M"
basepath = "cluster-results/sost-20210711-run-50M-g"

with open(os.path.join(basepath, "calculation.json"), 'r') as f:
    calcsettings = json.load(f)
    f.close()

batches = 10
sp = openmc.StatePoint(os.path.join(basepath, "statepoint.{:d}.h5".format(batches)))
tally = sp.get_tally(scores=['flux'])
#flux = tally.get_slice(scores=['flux'])

xwidth = len(tally.filters[0].mesh.x_grid) - 1
ywidth = len(tally.filters[0].mesh.y_grid) - 1
zwidth = len(tally.filters[0].mesh.z_grid) - 1
fluxdata = tally.get_reshaped_data()[:,0,0]
# fluxdata.shape = (xwidth, ywidth, zwidth)
fluxdata = fluxdata.reshape((zwidth, ywidth, xwidth))

plotcols = math.ceil(math.sqrt(len(fluxdata)))
plotrows = math.ceil(len(fluxdata) / plotcols)
fix, ax = plt.subplots(plotcols, plotrows, figsize=(15, 15))

# thanks to https://stackoverflow.com/questions/9455044/problems-with-zeros-in-matplotlib-colors-lognorm
my_cmap = copy.copy(matplotlib.cm.get_cmap('viridis')) # copy the default cmap
my_cmap.set_bad(my_cmap(0))

for i in range(len(fluxdata)):
    colidx = i % plotrows
    rowidx = i // plotcols
    print(rowidx, colidx)
    ax[colidx, rowidx].imshow(fluxdata[i], norm=LogNorm(vmin=0.000001, vmax=10), cmap=my_cmap)
    ax[colidx, rowidx].title.set_text("z={:d}, Max: {:.3e}".format(i, fluxdata[i].max()))

plt.savefig("soviet-storage-neutron-flux.png")    
plt.show()

i = 3
fig, ax = plt.subplots(figsize=(5, 5))
ax.imshow(fluxdata[i], norm=LogNorm(vmin=0.000001, vmax=10), cmap=my_cmap)
ax.title.set_text("z={:d}, Max: {:.3e}".format(i, fluxdata[i].max()))
plt.savefig("soviet-storage-neutron-flux-3m.png")    
plt.show()

i = 15
fig, ax = plt.subplots(figsize=(5, 5))
ax.imshow(fluxdata[i], norm=LogNorm(vmin=1e-13, vmax=1e-6), cmap=my_cmap)
ax.title.set_text("z={:d}, Max: {:.3e}".format(i, fluxdata[i].max()))
plt.savefig("soviet-storage-neutron-flux-3m.png")    
plt.show()

plotcols = math.ceil(math.sqrt(len(fluxdata)))
plotrows = math.ceil(len(fluxdata) / plotcols)
fix, ax = plt.subplots(plotcols, plotrows, figsize=(15, 15))
for i in range(len(fluxdata)):
    colidx = i % plotrows
    rowidx = i // plotcols
    print(rowidx, colidx)
    ax[colidx, rowidx].imshow(fluxdata[i], norm=LogNorm(vmin=1e-13, vmax=1e-6), cmap=my_cmap)
    ax[colidx, rowidx].title.set_text("z={:d}, Max: {:.3e}".format(i, fluxdata[i].max()))

plt.savefig("soviet-storage-neutron-flux.png")    
plt.show()
