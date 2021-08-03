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

import matplotlib.patches as patches
from PIL import Image

import helper

if len(sys.argv) == 2:
    basepath = sys.argv[1]
    cosmicbasepath = ""
else:
    basepath = "cluster-results/sost-20210730-fetter-survival-50M"
    cosmicbasepath = "cluster-results/sost-20210730-cosmic-200MeV-discard-20M"

pumass = 4000
cosmicneutronsperm2 = 120
days = 2

plotpath = os.path.join("plots", "sost")
if not os.path.exists(plotpath):
    os.mkdir(plotpath)

# thanks to https://stackoverflow.com/questions/9455044/problems-with-zeros-in-matplotlib-colors-lognorm
my_cmap = copy.copy(matplotlib.cm.get_cmap('viridis')) # copy the default cmap
# my_cmap.set_bad(my_cmap(0))
my_cmap.set_bad("black")
    
################################################################################
# Load weapon source statepoint file, flux and current
with open(os.path.join(basepath, "calculation.json"), 'r') as f:
    calcsettings = json.load(f)
    batches = calcsettings['batches']
    f.close()

sp = openmc.StatePoint(os.path.join(basepath, "statepoint.{:d}.h5".format(batches)))


puvec = helper.createwgpuvec()
age = calcsettings['weaponage']
puvec.createagedvector([age])
sourcesum = 0
for iso in puvec.pudf.index:
    sourcesum += puvec.puagedf.loc[(iso, age), 'sfneutrons'] * puvec.puagedf.loc[(iso, age), 'wo']
nps = sourcesum * pumass
print("The WPu weapon source ({:.2f} years old, {:.1f}g) emits {:.1f} neutrons/s".format(age, pumass, nps))

tally = sp.get_tally(name = 'flux, regular mesh')
(xwidth, ywidth, zwidth) = tally.find_filter(openmc.MeshFilter).mesh.dimension
print(xwidth, ywidth)

fluxdata = tally.get_reshaped_data()[:,0,0]
fluxdata = fluxdata.reshape((zwidth, ywidth, xwidth))
fluxdata *= nps

tally = sp.get_tally(name = 'current, regular mesh')
(xwidth, ywidth, zwidth) = tally.find_filter(openmc.MeshSurfaceFilter).mesh.dimension
meshid = tally.find_filter(openmc.MeshSurfaceFilter).mesh.id
meshstring = 'mesh {:d}'.format(meshid)
currentdf = tally.get_pandas_dataframe()
currentdf['mean'] = currentdf['mean'] * nps
currentdf['std. dev.'] = currentdf['std. dev.'] * nps

incurrentdf = currentdf[currentdf[(meshstring, 'surf')].isin(['x-min in', 'x-max in', 'y-min in', 'y-max in', 'z-min in', 'z-max in'])]
incurrentdf = incurrentdf.groupby([(meshstring, 'x'), (meshstring, 'y'), (meshstring, 'z')]).sum()
incurrentdf.reset_index(inplace = True)

currentdata = incurrentdf.sort_values([(meshstring, 'z'), (meshstring, 'y'), (meshstring, 'x')])['mean'].to_numpy()
currentdata = currentdata.reshape((zwidth, ywidth, xwidth))

shieldedincurrentdf = currentdf[currentdf[(meshstring, 'surf')].isin(['x-min in', 'x-max in', 'y-min in', 'y-max in', 'z-min in', 'z-max in'])]
shieldedincurrentdf = shieldedincurrentdf.groupby([(meshstring, 'x'), (meshstring, 'y'), (meshstring, 'z')]).sum()
shieldedincurrentdf.reset_index(inplace = True)

shieldedcurrentdata = shieldedincurrentdf.sort_values([(meshstring, 'z'), (meshstring, 'y'), (meshstring, 'x')])['mean'].to_numpy()
shieldedcurrentdata = currentdata.reshape((zwidth, ywidth, xwidth))


################################################################################
# Load cosmic ray statepoint, flux and current
with open(os.path.join(cosmicbasepath, "calculation.json"), 'r') as f:
    cosmiccalcsettings = json.load(f)
    cosmicbatches = cosmiccalcsettings['batches']
    f.close()

cosmicsp = openmc.StatePoint(os.path.join(cosmicbasepath, "statepoint.{:d}.h5".format(cosmicbatches)))

if cosmicbasepath != "":
    cosmicsettings = openmc.Settings.from_xml(os.path.join(cosmicbasepath, "settings.xml"))
    parameters = cosmicsettings.source[0].parameters
    cosmicraywidth = float(parameters.split("subboxLength ")[-1].split(" ")[0])
    cosmicneutrons = cosmicneutronsperm2 * cosmicraywidth ** 2
else:
    cosmicneutrons = cosmicneutronsperm2

cosmictally = cosmicsp.get_tally(name = 'flux, regular mesh')
(cosmicxwidth, cosmicywidth, cosmiczwidth) = cosmictally.find_filter(openmc.MeshFilter).mesh.dimension
print(cosmicxwidth, cosmicywidth)

cosmicfluxdata = cosmictally.get_reshaped_data()[:,0,0]
cosmicfluxdata = cosmicfluxdata.reshape((cosmiczwidth, cosmicywidth, cosmicxwidth))
cosmicfluxdata *= cosmicneutrons

cosmictally = cosmicsp.get_tally(name = 'current, regular mesh')
cosmicmeshid = cosmictally.find_filter(openmc.MeshSurfaceFilter).mesh.id
cosmicmeshstring = 'mesh {:d}'.format(cosmicmeshid)
(cosmicxwidth, cosmicywidth, cosmiczwidth) = cosmictally.find_filter(openmc.MeshSurfaceFilter).mesh.dimension
print(cosmicxwidth, cosmicywidth)

cosmiccurrentdf = cosmictally.get_pandas_dataframe()
cosmiccurrentdf['mean'] = cosmiccurrentdf['mean'] * cosmicneutrons
cosmiccurrentdf['std. dev.'] = cosmiccurrentdf['std. dev.'] * cosmicneutrons
cosmicincurrentdf = cosmiccurrentdf[cosmiccurrentdf[(cosmicmeshstring, 'surf')].isin(['x-min in', 'x-max in', 'y-min in', 'y-max in', 'z-min in', 'z-max in'])]
cosmicincurrentdf = cosmicincurrentdf.groupby([(cosmicmeshstring, 'x'), (cosmicmeshstring, 'y'), (cosmicmeshstring, 'z')]).sum()
cosmicincurrentdf.reset_index(inplace = True)

cosmiccurrentdata = cosmicincurrentdf.sort_values([(cosmicmeshstring, 'z'), (cosmicmeshstring, 'y'), (cosmicmeshstring, 'x')])['mean'].to_numpy()
cosmiccurrentdata = cosmiccurrentdata.reshape((cosmiczwidth, cosmicywidth, cosmicxwidth))

shieldedcosmicincurrentdf = cosmiccurrentdf[cosmiccurrentdf[(cosmicmeshstring, 'surf')].isin(['x-min in', 'x-max in', 'y-min in', 'y-max in'])]
shieldedcosmicincurrentdf = shieldedcosmicincurrentdf.groupby([(cosmicmeshstring, 'x'), (cosmicmeshstring, 'y'), (cosmicmeshstring, 'z')]).sum()
shieldedcosmicincurrentdf.reset_index(inplace = True)

shieldedcosmiccurrentdata = shieldedcosmicincurrentdf.sort_values([(cosmicmeshstring, 'z'), (cosmicmeshstring, 'y'), (cosmicmeshstring, 'x')])['mean'].to_numpy()
shieldedcosmiccurrentdata = shieldedcosmiccurrentdata.reshape((cosmiczwidth, cosmicywidth, cosmicxwidth))

tempdf = cosmictally.get_pandas_dataframe()
tempdf = tempdf[tempdf[(cosmicmeshstring, 'surf')].isin(['z-max in'])]
currenttowardsearth = tempdf[tempdf[(cosmicmeshstring, 'z')] == 15]['mean'].sum()
print("total current towards earth z-min out, at z = 15m: {:f}".format(currenttowardsearth))
print("total cosmic current (all in), at z = 14m: {:f}".format(cosmiccurrentdata[-2].sum() / cosmicneutrons))

################################################################################
# Analysis
groundfloorlevel = sp.summary.geometry.get_all_surfaces()[2020].z0
groundfloormeshid = int(math.ceil(groundfloorlevel / 100))
# groundfloormeshid = 4
print("The ground is at z={:.2f}cm".format(groundfloorlevel))
print("A detector on the ground will be in mesh cell with z-index {:d}".format(groundfloormeshid))

################################################################################
# Plot measurement times for 3 std dev (shielded) one level
m = 3

vmax = 24 * 3600 * days

for z in range(zwidth):
    fig, ax = plt.subplots(nrows = 1, ncols = 1, squeeze = False, figsize = (5, 6))
    data = m ** 2 * shieldedcosmiccurrentdata[z] / (shieldedcurrentdata[z] ** 2)
    data[data == np.inf] = 0
    im = ax[0, 0].imshow(data, norm=LogNorm(vmax = vmax), cmap=my_cmap)
    ax[0, 0].title.set_text("z={:d}, Max: {:.3e}".format(z, data.max()))
    # fig.colorbar(im, ax = ax[0, 0])
    # axis need to be in list
    # https://stackoverflow.com/questions/37280557/can-i-place-a-vertical-colorbar-to-the-left-of-the-plot-in-matplotlib    
    fig.colorbar(im, ax = [ax[0, 0]], location = 'bottom')

    plt.xlim(60, 140)
    plt.ylim(60, 140)
    plt.savefig(os.path.join(plotpath, "soviet-storage-z-{}-shielded-measurement-time-{}stddev.png".format(z, m)))
    plt.savefig(os.path.join(plotpath, "soviet-storage-z-{}-shielded-measurement-time-{}stddev.pdf".format(z, m)))
    # plt.show()
    plt.close()

################################################################################
# Plot shielded current from weapon, one level
# vmin = 0.0001
# vmax = 10 ** (math.ceil(math.log(currentdata.max(), 10)) - 1)

for z in range(zwidth):
    fig, ax = plt.subplots(nrows = 1, ncols = 1, squeeze = False, figsize = (5, 6))
    data = shieldedcurrentdata[z]
    data[data == np.inf] = 0
    im = ax[0, 0].imshow(data, norm=LogNorm(), cmap=my_cmap)
    ax[0, 0].title.set_text("z={:d}, Max: {:.3e}".format(z, data.max()))
    # fig.colorbar(im, ax = ax[0, 0])
    # axis need to be in list
    # https://stackoverflow.com/questions/37280557/can-i-place-a-vertical-colorbar-to-the-left-of-the-plot-in-matplotlib    
    fig.colorbar(im, ax = [ax[0, 0]], location = 'bottom')

    plt.xlim(60, 140)
    plt.ylim(60, 140)
    plt.savefig(os.path.join(plotpath, "soviet-storage-z-{}-shielded-current.png".format(z)))
    plt.savefig(os.path.join(plotpath, "soviet-storage-z-{}-shielded-current.pdf".format(z)))
    # plt.show()
    plt.close()

################################################################################
# Plot shielded current from cosmic rays, one level
# vmin = 0.0001
# vmax = 10 ** (math.ceil(math.log(currentdata.max(), 10)) - 1)

for z in range(zwidth):
    fig, ax = plt.subplots(nrows = 1, ncols = 1, squeeze = False, figsize = (5, 6))
    data = shieldedcosmiccurrentdata[z]
    data[data == np.inf] = 0
    im = ax[0, 0].imshow(data, norm=LogNorm(), cmap=my_cmap)
    ax[0, 0].title.set_text("z={:d}, Max: {:.3e}".format(z, data.max()))
    # fig.colorbar(im, ax = ax[0, 0])
    # axis need to be in list
    # https://stackoverflow.com/questions/37280557/can-i-place-a-vertical-colorbar-to-the-left-of-the-plot-in-matplotlib    
    fig.colorbar(im, ax = [ax[0, 0]], location = 'bottom')

    plt.xlim(60, 140)
    plt.ylim(60, 140)
    plt.savefig(os.path.join(plotpath, "soviet-storage-z-{}-shielded-cosmic-current.png".format(z)))
    plt.savefig(os.path.join(plotpath, "soviet-storage-z-{}-shielded-cosmic-current.pdf".format(z)))
    # plt.show()
    plt.close()

################################################################################
# Plot flux of weapon source
vmin = 1e-2
vmax = 10 ** math.ceil(math.log(nps, 10))

plotcols = math.ceil(math.sqrt(len(fluxdata)))
plotrows = math.ceil(len(fluxdata) / plotcols)
fig, ax = plt.subplots(nrows = plotrows, ncols = plotcols, squeeze = False,
                       figsize=(20, 20), sharex = True, sharey = True)

for i in range(len(fluxdata)):
    colidx = i % plotrows
    rowidx = i // plotcols
    print(i, rowidx, colidx)
    if i < len(fluxdata):
        im = ax[rowidx, colidx].imshow(fluxdata[i], norm=LogNorm(vmin=vmin, vmax=vmax), cmap=my_cmap)
        ax[rowidx, colidx].title.set_text("z={:d}, Max: {:.3e}".format(i, fluxdata[i].max()))
    else:
        fig.delaxes(ax[rowidx, colidx])

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(im, cax=cbar_ax)

plt.savefig(os.path.join(plotpath, "soviet-storage-neutron-flux.png"))
#plt.show()
plt.close()

################################################################################
# Plot measurement times for 1 std dev
m = 3

plotcols = math.ceil(math.sqrt(len(fluxdata)))
plotrows = math.ceil(len(fluxdata) / plotcols)
fig, ax = plt.subplots(nrows = plotrows, ncols = plotcols, squeeze = False,
                       figsize=(20, 20), sharex = True, sharey = True)

vmax = 24 * 3600 * days

for i in range(len(fluxdata)):
    colidx = i % plotrows
    rowidx = i // plotcols
    print(i, rowidx, colidx)
    data = m ** 2 * cosmiccurrentdata[i] / (currentdata[i] ** 2)
    data[data == np.inf] = 0
    if i < len(fluxdata):
        # im = ax[rowidx, colidx].imshow(fluxdata[i], norm=LogNorm(vmin=vmin, vmax=vmax), cmap=my_cmap)
        im = ax[rowidx, colidx].imshow(data, norm=LogNorm(vmax = vmax), cmap=my_cmap)
        ax[rowidx, colidx].title.set_text("z={:d}, Max: {:.3e}".format(i, data.max()))
        fig.colorbar(im, ax = ax[rowidx, colidx])
    else:
        fig.delaxes(ax[rowidx, colidx])

plt.xlim(60, 140)
plt.ylim(60, 140)
plt.savefig(os.path.join(plotpath, "soviet-storage-measurement-time-{}stddev.png".format(m)))
plt.show()
plt.close()

################################################################################
# Plot measurement times for 1 std dev (shielded)
m = 3

plotcols = math.ceil(math.sqrt(len(fluxdata)))
plotrows = math.ceil(len(fluxdata) / plotcols)
fig, ax = plt.subplots(nrows = plotrows, ncols = plotcols, squeeze = False,
                       figsize=(20, 20), sharex = True, sharey = True)

vmax = 24 * 3600 * days

for i in range(len(fluxdata)):
    colidx = i % plotrows
    rowidx = i // plotcols
    print(i, rowidx, colidx)
    data = m ** 2 * shieldedcosmiccurrentdata[i] / (shieldedcurrentdata[i] ** 2)
    data[data == np.inf] = 0
    if i < len(fluxdata):
        # im = ax[rowidx, colidx].imshow(fluxdata[i], norm=LogNorm(vmin=vmin, vmax=vmax), cmap=my_cmap)
        im = ax[rowidx, colidx].imshow(data, norm=LogNorm(vmax = vmax), cmap=my_cmap)
        ax[rowidx, colidx].title.set_text("z={:d}, Max: {:.3e}".format(i, data.max()))
        fig.colorbar(im, ax = ax[rowidx, colidx])
    else:
        fig.delaxes(ax[rowidx, colidx])

plt.xlim(60, 140)
plt.ylim(60, 140)
plt.savefig(os.path.join(plotpath, "soviet-storage-shielded-measurement-time-{}stddev.png".format(m)))
plt.show()
plt.close()


################################################################################
# Plot current at ground floor, compare weapon source and cosmic rays

vmin = 0.1
vmax = 10 ** (math.ceil(math.log(currentdata.max(), 10)) - 1)

fig, ax = plt.subplots(1, 2, squeeze = False, sharey = True, figsize=(20, 10))
im = ax[0, 0].imshow(currentdata[groundfloormeshid], norm=LogNorm(vmin, vmax), cmap=my_cmap)
ax[0, 1].imshow(cosmiccurrentdata[groundfloormeshid], norm=LogNorm(vmin, vmax), cmap=my_cmap)
plt.show()

x = math.floor(xwidth / 2)
cosmicx = math.floor(cosmicxwidth / 2)
y = 17
cosmicy = 17

def currentcompare(x, y, z):
    xoffset = (cosmicxwidth - xwidth) // 2
    wx = x + xoffset
    print(wx, x)
    print("Weapon current at ({:d}, {:d}, {:d}): {:e}".format(wx, y, z, currentdata[z, y, wx]))
    print("Cosmic ray current at ({:d}, {:d}, {:d}): {:e}".format(x, y, z, cosmiccurrentdata[z, y, x]))

    boxdf = currentdf[(currentdf[(meshstring, 'x')] == wx + 1) &
                      (currentdf[(meshstring, 'y')] == y + 1) &
                      (currentdf[(meshstring, 'z')] == z + 1)]
    cosmicboxdf = cosmiccurrentdf[(cosmiccurrentdf[(meshstring, 'x')] == x + 1) &
                                  (cosmiccurrentdf[(meshstring, 'y')] == y + 1) &
                                  (cosmiccurrentdf[(meshstring, 'z')] == z + 1)]
    print(boxdf)
    print(cosmicboxdf)

currentcompare(x, y, groundfloormeshid)

