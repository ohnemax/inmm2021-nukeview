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
from mpl_toolkits.axes_grid1 import make_axes_locatable

import helper

bustbp = "cluster-results/bust-20210801-fetter-survival-0.05-200M"
bustcosmicbp = "cluster-results/bust-20210801-cosmic-200MeV-discard-50M"
sostbp = "cluster-results/sost-20210801-fetter-survival-0.05-200M"
sostcosmicbp = "cluster-results/sost-20210801-cosmic-200MeV-discard-50M"

pumass = 4000
cosmicneutronsperm2 = 120
days = 1

# thanks to https://stackoverflow.com/questions/9455044/problems-with-zeros-in-matplotlib-colors-lognorm
my_cmap = copy.copy(matplotlib.cm.get_cmap('viridis')) # copy the default cmap
# my_cmap.set_bad(my_cmap(0))
my_cmap.set_bad("black")

puvec = helper.createwgpuvec()
age = 0
puvec.createagedvector([age])
sourcesum = 0
for iso in puvec.pudf.index:
    sourcesum += puvec.puagedf.loc[(iso, age), 'sfneutrons'] * puvec.puagedf.loc[(iso, age), 'wo']
nps = sourcesum * pumass
print("The WPu weapon source ({:.2f} years old, {:.1f}g) emits {:.1f} neutrons/s".format(age, pumass, nps))

################################################################################
# Büchel: Load weapon source statepoint file, flux and current
print("Loading Büchel, weapon source")
with open(os.path.join(bustbp, "calculation.json"), 'r') as f:
    calcsettings = json.load(f)
    batches = calcsettings['batches']
    f.close()
sp = openmc.StatePoint(os.path.join(bustbp, "statepoint.{:d}.h5".format(batches)))

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

bustcurrentdata = incurrentdf.sort_values([(meshstring, 'z'), (meshstring, 'y'), (meshstring, 'x')])['mean'].to_numpy()
bustcurrentdata = bustcurrentdata.reshape((zwidth, ywidth, xwidth))

shieldedincurrentdf = currentdf[currentdf[(meshstring, 'surf')].isin(['x-min in', 'x-max in', 'y-min in', 'y-max in', 'z-min in', 'z-max in'])]
shieldedincurrentdf = shieldedincurrentdf.groupby([(meshstring, 'x'), (meshstring, 'y'), (meshstring, 'z')]).sum()
shieldedincurrentdf.reset_index(inplace = True)

bustshieldedcurrentdata = shieldedincurrentdf.sort_values([(meshstring, 'z'), (meshstring, 'y'), (meshstring, 'x')])['mean'].to_numpy()
bustshieldedcurrentdata = bustshieldedcurrentdata.reshape((zwidth, ywidth, xwidth))

################################################################################
# Büchel: Load cosmic ray statepoint, flux and current
print("Loading Büchel, cosmic source")
with open(os.path.join(bustcosmicbp, "calculation.json"), 'r') as f:
    cosmiccalcsettings = json.load(f)
    cosmicbatches = cosmiccalcsettings['batches']
    f.close()

cosmicsp = openmc.StatePoint(os.path.join(bustcosmicbp, "statepoint.{:d}.h5".format(cosmicbatches)))

cosmicsettings = openmc.Settings.from_xml(os.path.join(bustcosmicbp, "settings.xml"))
parameters = cosmicsettings.source[0].parameters
bustcosmicraywidth = float(parameters.split("subboxLength ")[-1].split(" ")[0])
bustcosmicneutrons = cosmicneutronsperm2 * bustcosmicraywidth ** 2

cosmictally = cosmicsp.get_tally(name = 'current, regular mesh')
cosmicmeshid = cosmictally.find_filter(openmc.MeshSurfaceFilter).mesh.id
cosmicmeshstring = 'mesh {:d}'.format(cosmicmeshid)
(cosmicxwidth, cosmicywidth, cosmiczwidth) = cosmictally.find_filter(openmc.MeshSurfaceFilter).mesh.dimension
print(cosmicxwidth, cosmicywidth)

cosmiccurrentdf = cosmictally.get_pandas_dataframe()
cosmiccurrentdf['mean'] = cosmiccurrentdf['mean'] * bustcosmicneutrons
cosmiccurrentdf['std. dev.'] = cosmiccurrentdf['std. dev.'] * bustcosmicneutrons
cosmicincurrentdf = cosmiccurrentdf[cosmiccurrentdf[(cosmicmeshstring, 'surf')].isin(['x-min in', 'x-max in', 'y-min in', 'y-max in', 'z-min in', 'z-max in'])]
cosmicincurrentdf = cosmicincurrentdf.groupby([(cosmicmeshstring, 'x'), (cosmicmeshstring, 'y'), (cosmicmeshstring, 'z')]).sum()
cosmicincurrentdf.reset_index(inplace = True)

bustcosmiccurrentdata = cosmicincurrentdf.sort_values([(cosmicmeshstring, 'z'), (cosmicmeshstring, 'y'), (cosmicmeshstring, 'x')])['mean'].to_numpy()
bustcosmiccurrentdata = bustcosmiccurrentdata.reshape((cosmiczwidth, cosmicywidth, cosmicxwidth))

shieldedcosmicincurrentdf = cosmiccurrentdf[cosmiccurrentdf[(cosmicmeshstring, 'surf')].isin(['x-min in', 'x-max in', 'y-min in', 'y-max in'])]
shieldedcosmicincurrentdf = shieldedcosmicincurrentdf.groupby([(cosmicmeshstring, 'x'), (cosmicmeshstring, 'y'), (cosmicmeshstring, 'z')]).sum()
shieldedcosmicincurrentdf.reset_index(inplace = True)

bustshieldedcosmiccurrentdata = shieldedcosmicincurrentdf.sort_values([(cosmicmeshstring, 'z'), (cosmicmeshstring, 'y'), (cosmicmeshstring, 'x')])['mean'].to_numpy()
bustshieldedcosmiccurrentdata = bustshieldedcosmiccurrentdata.reshape((cosmiczwidth, cosmicywidth, cosmicxwidth))

################################################################################
# Soviet Storage: Load weapon source statepoint file, flux and current
print("Loading Soviet Storage, weapon source")
with open(os.path.join(sostbp, "calculation.json"), 'r') as f:
    calcsettings = json.load(f)
    batches = calcsettings['batches']
    f.close()
sp = openmc.StatePoint(os.path.join(sostbp, "statepoint.{:d}.h5".format(batches)))

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

sostcurrentdata = incurrentdf.sort_values([(meshstring, 'z'), (meshstring, 'y'), (meshstring, 'x')])['mean'].to_numpy()
sostcurrentdata = sostcurrentdata.reshape((zwidth, ywidth, xwidth))

shieldedincurrentdf = currentdf[currentdf[(meshstring, 'surf')].isin(['x-min in', 'x-max in', 'y-min in', 'y-max in', 'z-min in', 'z-max in'])]
shieldedincurrentdf = shieldedincurrentdf.groupby([(meshstring, 'x'), (meshstring, 'y'), (meshstring, 'z')]).sum()
shieldedincurrentdf.reset_index(inplace = True)

sostshieldedcurrentdata = shieldedincurrentdf.sort_values([(meshstring, 'z'), (meshstring, 'y'), (meshstring, 'x')])['mean'].to_numpy()
sostshieldedcurrentdata = sostshieldedcurrentdata.reshape((zwidth, ywidth, xwidth))

################################################################################
# Soviet Storage: Load cosmic ray statepoint, flux and current
print("Loading Soviet Storage, cosmic source")
with open(os.path.join(sostcosmicbp, "calculation.json"), 'r') as f:
    cosmiccalcsettings = json.load(f)
    cosmicbatches = cosmiccalcsettings['batches']
    f.close()

cosmicsp = openmc.StatePoint(os.path.join(sostcosmicbp, "statepoint.{:d}.h5".format(cosmicbatches)))

cosmicsettings = openmc.Settings.from_xml(os.path.join(sostcosmicbp, "settings.xml"))
parameters = cosmicsettings.source[0].parameters
sostcosmicraywidth = float(parameters.split("subboxLength ")[-1].split(" ")[0])
sostcosmicneutrons = cosmicneutronsperm2 * sostcosmicraywidth ** 2

cosmictally = cosmicsp.get_tally(name = 'current, regular mesh')
cosmicmeshid = cosmictally.find_filter(openmc.MeshSurfaceFilter).mesh.id
cosmicmeshstring = 'mesh {:d}'.format(cosmicmeshid)
(cosmicxwidth, cosmicywidth, cosmiczwidth) = cosmictally.find_filter(openmc.MeshSurfaceFilter).mesh.dimension
print(cosmicxwidth, cosmicywidth)

cosmiccurrentdf = cosmictally.get_pandas_dataframe()
cosmiccurrentdf['mean'] = cosmiccurrentdf['mean'] * sostcosmicneutrons
cosmiccurrentdf['std. dev.'] = cosmiccurrentdf['std. dev.'] * sostcosmicneutrons
cosmicincurrentdf = cosmiccurrentdf[cosmiccurrentdf[(cosmicmeshstring, 'surf')].isin(['x-min in', 'x-max in', 'y-min in', 'y-max in', 'z-min in', 'z-max in'])]
cosmicincurrentdf = cosmicincurrentdf.groupby([(cosmicmeshstring, 'x'), (cosmicmeshstring, 'y'), (cosmicmeshstring, 'z')]).sum()
cosmicincurrentdf.reset_index(inplace = True)

sostcosmiccurrentdata = cosmicincurrentdf.sort_values([(cosmicmeshstring, 'z'), (cosmicmeshstring, 'y'), (cosmicmeshstring, 'x')])['mean'].to_numpy()
sostcosmiccurrentdata = sostcosmiccurrentdata.reshape((cosmiczwidth, cosmicywidth, cosmicxwidth))

shieldedcosmicincurrentdf = cosmiccurrentdf[cosmiccurrentdf[(cosmicmeshstring, 'surf')].isin(['x-min in', 'x-max in', 'y-min in', 'y-max in'])]
shieldedcosmicincurrentdf = shieldedcosmicincurrentdf.groupby([(cosmicmeshstring, 'x'), (cosmicmeshstring, 'y'), (cosmicmeshstring, 'z')]).sum()
shieldedcosmicincurrentdf.reset_index(inplace = True)

sostshieldedcosmiccurrentdata = shieldedcosmicincurrentdf.sort_values([(cosmicmeshstring, 'z'), (cosmicmeshstring, 'y'), (cosmicmeshstring, 'x')])['mean'].to_numpy()
sostshieldedcosmiccurrentdata = sostshieldedcosmiccurrentdata.reshape((cosmiczwidth, cosmicywidth, cosmicxwidth))

################################################################################
# Combined plot
m = 3
bustplotz = 3
sostplotz1 = 2
sostplotz2 = 7
pad = 5
colorbarsize = "8%"

tmin = 1
tmax = 24 * 3600

factor = 0.8
fig, ax = plt.subplots(nrows = 3, ncols = 3, squeeze = False,
                       figsize = (12 * factor, 3 * 3 * factor)
)

vmin = 0.0001
vmax = 10 ** (math.ceil(math.log(bustshieldedcosmiccurrentdata.max(), 10)))

data = bustshieldedcurrentdata[bustplotz]
im = ax[0, 0].imshow(data,
                     norm=LogNorm(vmin = vmin, vmax = vmax),
                     cmap=my_cmap,
                     interpolation = 'none')
divider = make_axes_locatable(ax[0, 0])
cax1 = divider.append_axes("right", size=colorbarsize, pad=0.05)
fig.colorbar(im, cax = cax1,
             label = "Inward Current [neutron/s]")

data = bustshieldedcosmiccurrentdata[bustplotz]
im = ax[0, 1].imshow(data,
                     norm=LogNorm(vmin = vmin, vmax = vmax),
                     cmap=my_cmap,
                     interpolation = 'none')
divider = make_axes_locatable(ax[0, 1])
cax2 = divider.append_axes("right", size=colorbarsize, pad=0.05)
fig.colorbar(im, cax = cax2,
             label = "Inward Current [neutron/s]")

data = m ** 2 * bustshieldedcosmiccurrentdata[bustplotz] / (bustshieldedcurrentdata[bustplotz] ** 2)
data[data == np.inf] = 0
im = ax[0, 2].imshow(data,
                     norm=LogNorm(vmin = tmin, vmax = tmax),
                     cmap=my_cmap,
                     interpolation = 'none')
divider = make_axes_locatable(ax[0, 2])
cax3 = divider.append_axes("right", size=colorbarsize, pad=0.05)
fig.colorbar(im, cax = cax3,
             label = "Measurement Time [s]")

ax[0, 0].set_xlim(70, 130)
ax[0, 0].set_ylim(70, 130)
ax[0, 1].set_xlim(70, 130)
ax[0, 1].set_ylim(70, 130)
ax[0, 2].set_xlim(70, 130)
ax[0, 2].set_ylim(70, 130)

ax[0, 0].set_ylabel("y-axis [m]")

ax[0, 0].annotate("Büchel PSA", xy=(0, 0.5), xytext=(-ax[0, 0].yaxis.labelpad - pad, 0),
                  xycoords=ax[0, 0].yaxis.label, textcoords='offset points',
                  size='large', ha='right', va='center', rotation = 90)

vmin = 0.0001
vmax = 10 ** (math.ceil(math.log(sostshieldedcosmiccurrentdata.max(), 10)))

data = sostshieldedcurrentdata[sostplotz1]
im = ax[1, 0].imshow(data,
                     norm=LogNorm(vmin = vmin, vmax = vmax),
                     cmap=my_cmap,
                     interpolation = 'none')
divider = make_axes_locatable(ax[1, 0])
cax1 = divider.append_axes("right", size=colorbarsize, pad=0.05)
fig.colorbar(im, cax = cax1,
             label = "Inward Current [neutron/s]")

data = sostshieldedcosmiccurrentdata[sostplotz1]
im = ax[1, 1].imshow(data,
                     norm=LogNorm(vmin = vmin, vmax = vmax),
                     cmap=my_cmap,
                     interpolation = 'none')
divider = make_axes_locatable(ax[1, 1])
cax2 = divider.append_axes("right", size=colorbarsize, pad=0.05)
fig.colorbar(im, cax = cax2,
             label = "Inward Current [neutron/s]")

data = m ** 2 * sostshieldedcosmiccurrentdata[sostplotz1] / (sostshieldedcurrentdata[sostplotz1] ** 2)
data[data == np.inf] = 0
im = ax[1, 2].imshow(data,
                     norm=LogNorm(vmin = tmin, vmax = tmax),
                     cmap=my_cmap,
                     interpolation = 'none')
divider = make_axes_locatable(ax[1, 2])
cax3 = divider.append_axes("right", size=colorbarsize, pad=0.05)
fig.colorbar(im, cax = cax3,
             label = "Measurement Time [s]")

ax[1, 0].set_xlim(60, 140)
ax[1, 0].set_ylim(60, 140)
ax[1, 1].set_xlim(60, 140)
ax[1, 1].set_ylim(60, 140)
ax[1, 2].set_xlim(60, 140)
ax[1, 2].set_ylim(60, 140)

ax[1, 0].set_ylabel("y-axis [m]")

ax[1, 0].annotate("Soviet Storage, {}m to {}m".format(sostplotz1, sostplotz1 + 1), xy=(0, 0.5), xytext=(-ax[1, 0].yaxis.labelpad - pad, 0),
                  xycoords=ax[1, 0].yaxis.label, textcoords='offset points',
                  size='large', ha='right', va='center', rotation = 90)

vmin = 0.0001
vmax = 10 ** (math.ceil(math.log(sostshieldedcosmiccurrentdata.max(), 10)))

data = sostshieldedcurrentdata[sostplotz2]
im = ax[2, 0].imshow(data,
                     norm=LogNorm(vmin = vmin, vmax = vmax),
                     cmap=my_cmap,
                     interpolation = 'none')
divider = make_axes_locatable(ax[2, 0])
cax1 = divider.append_axes("right", size=colorbarsize, pad=0.05)
fig.colorbar(im, cax = cax1,
             label = "Inward Current [neutron/s]")

data = sostshieldedcosmiccurrentdata[sostplotz2]
im = ax[2, 1].imshow(data,
                     norm=LogNorm(vmin = vmin, vmax = vmax),
                     cmap=my_cmap,
                     interpolation = 'none')
divider = make_axes_locatable(ax[2, 1])
cax2 = divider.append_axes("right", size=colorbarsize, pad=0.05)
fig.colorbar(im, cax = cax2,
             label = "Inward Current [neutron/s]")

data = m ** 2 * sostshieldedcosmiccurrentdata[sostplotz2] / (sostshieldedcurrentdata[sostplotz2] ** 2)
data[data == np.inf] = 0
im = ax[2, 2].imshow(data,
                     norm=LogNorm(vmin = tmin, vmax = tmax),
                     cmap=my_cmap,
                     interpolation = 'none')
divider = make_axes_locatable(ax[2, 2])
cax3 = divider.append_axes("right", size=colorbarsize, pad=0.05)
fig.colorbar(im, cax = cax3,
             label = "Measurement Time [s]")

ax[2, 0].set_xlim(60, 140)
ax[2, 0].set_ylim(60, 140)
ax[2, 1].set_xlim(60, 140)
ax[2, 1].set_ylim(60, 140)
ax[2, 2].set_xlim(60, 140)
ax[2, 2].set_ylim(60, 140)

ax[2, 0].set_ylabel("y-axis [m]")
ax[2, 0].set_xlabel("x-axis [m]")
ax[2, 1].set_xlabel("x-axis [m]")
ax[2, 2].set_xlabel("x-axis [m]")

ax[2, 0].annotate("Soviet Storage, {}m to {}m".format(sostplotz2, sostplotz2 + 1), xy=(0, 0.5), xytext=(-ax[2, 0].yaxis.labelpad - pad, 0),
                  xycoords=ax[2, 0].yaxis.label, textcoords='offset points',
                  size='large', ha='right', va='center', rotation = 90)

ax[0, 0].annotate("Source: Weapon", xy=(0.5, 1), xytext=(0, pad),
                  xycoords='axes fraction', textcoords='offset points',
                  size='large', ha='center', va='baseline')
ax[0, 1].annotate("Source: Cosmic", xy=(0.5, 1), xytext=(0, pad),
                  xycoords='axes fraction', textcoords='offset points',
                  size='large', ha='center', va='baseline')
ax[0, 2].annotate("Measurement Time", xy=(0.5, 1), xytext=(0, pad),
                  xycoords='axes fraction', textcoords='offset points',
                  size='large', ha='center', va='baseline')

fig.tight_layout()

plt.savefig(os.path.join("plots-for-paper", "buechel-soviet-combined-m{}.pdf".format(m)))
plt.savefig(os.path.join("plots-for-paper", "buechel-soviet-combined-m{}.png".format(m)))

plt.show()

################################################################################
# Combined plot Büchel for Presentation
m = 3
bustplotz = 3
pad = 0.5
titlepad = 5
colorbarsize = "8%"

tmin = 1
tmax = 24 * 3600

factor = 0.8
fig, ax = plt.subplots(nrows = 1, ncols = 3, squeeze = False,
                       figsize = (12 * factor, 5 * factor),
                       sharey = True
)

vmin = 0.0001
vmax = 10 ** (math.ceil(math.log(bustshieldedcosmiccurrentdata.max(), 10)))

data = bustshieldedcurrentdata[bustplotz]
im = ax[0, 0].imshow(data,
                     norm=LogNorm(vmin = vmin, vmax = vmax),
                     cmap=my_cmap,
                     interpolation = 'none')
divider = make_axes_locatable(ax[0, 0])
cax1 = divider.append_axes("bottom", size=colorbarsize, pad=pad)
fig.colorbar(im, cax = cax1,
             label = "Current into 1 m$^3$ [neutron/s]",
             orientation="horizontal")

data = bustshieldedcosmiccurrentdata[bustplotz]
im = ax[0, 1].imshow(data,
                     norm=LogNorm(vmin = vmin, vmax = vmax),
                     cmap=my_cmap,
                     interpolation = 'none')
divider = make_axes_locatable(ax[0, 1])
cax2 = divider.append_axes("bottom", size=colorbarsize, pad=pad)
fig.colorbar(im, cax = cax2,
             label = "Current into 1 m$^3$ [neutron/s]",
             orientation="horizontal")

data = m ** 2 * bustshieldedcosmiccurrentdata[bustplotz] / (bustshieldedcurrentdata[bustplotz] ** 2)
data[data == np.inf] = 0
im = ax[0, 2].imshow(data,
                     norm=LogNorm(vmin = tmin, vmax = tmax),
                     cmap=my_cmap,
                     interpolation = 'none')
divider = make_axes_locatable(ax[0, 2])
cax3 = divider.append_axes("bottom", size=colorbarsize, pad=pad)
fig.colorbar(im, cax = cax3,
             label = "Measurement Time [s]",
             orientation="horizontal")

ax[0, 0].set_xlim(70, 130)
ax[0, 0].set_ylim(70, 130)
ax[0, 1].set_xlim(70, 130)
ax[0, 1].set_ylim(70, 130)
ax[0, 2].set_xlim(70, 130)
ax[0, 2].set_ylim(70, 130)

ax[0, 0].set_ylabel("[meters]")
ax[0, 0].set_xlabel("[meters]")
ax[0, 1].set_xlabel("[meters]")
ax[0, 2].set_xlabel("[meters]")

ax[0, 0].annotate("Source: Weapon", xy=(0.5, 1), xytext=(0, titlepad),
                  xycoords='axes fraction', textcoords='offset points',
                  size='large', ha='center', va='baseline')
ax[0, 1].annotate("Source: Cosmic", xy=(0.5, 1), xytext=(0, titlepad),
                  xycoords='axes fraction', textcoords='offset points',
                  size='large', ha='center', va='baseline')
ax[0, 2].annotate("Measurement Time", xy=(0.5, 1), xytext=(0, titlepad),
                  xycoords='axes fraction', textcoords='offset points',
                  size='large', ha='center', va='baseline')

#plt.tight_layout()
plt.savefig(os.path.join("plots-for-presentation", "buechel-m{}.pdf".format(m)), transparent=True)
plt.savefig(os.path.join("plots-for-presentation", "buechel-m{}.png".format(m)), transparent=True)

plt.show()

################################################################################
# Combined plot Soviet Storage for Presentation
m = 3
sostplotz1 = 2
sostplotz2 = 7
pad = 0.5
titlepad = 5
colorbarsize = "8%"

tmin = 1
tmax = 24 * 3600

factor = 0.8
fig, ax = plt.subplots(nrows = 1, ncols = 3, squeeze = False,
                       figsize = (12 * factor, 5 * factor),
                       sharey = True
)

vmin = 0.0001
vmax = 10 ** (math.ceil(math.log(sostshieldedcosmiccurrentdata.max(), 10)))

data = sostshieldedcurrentdata[sostplotz1]
im = ax[0, 0].imshow(data,
                     norm=LogNorm(vmin = vmin, vmax = vmax),
                     cmap=my_cmap,
                     interpolation = 'none')
divider = make_axes_locatable(ax[0, 0])
cax1 = divider.append_axes("bottom", size=colorbarsize, pad=pad)
fig.colorbar(im, cax = cax1,
             label = "Current into 1 m$^3$ [neutron/s]",
             orientation="horizontal")

data = sostshieldedcosmiccurrentdata[sostplotz1]
im = ax[0, 1].imshow(data,
                     norm=LogNorm(vmin = vmin, vmax = vmax),
                     cmap=my_cmap,
                     interpolation = 'none')
divider = make_axes_locatable(ax[0, 1])
cax2 = divider.append_axes("bottom", size=colorbarsize, pad=pad)
fig.colorbar(im, cax = cax2,
             label = "Current into 1 m$^3$ [neutron/s]",
             orientation="horizontal")

data = m ** 2 * sostshieldedcosmiccurrentdata[sostplotz1] / (sostshieldedcurrentdata[sostplotz1] ** 2)
data[data == np.inf] = 0
im = ax[0, 2].imshow(data,
                     norm=LogNorm(vmin = tmin, vmax = tmax),
                     cmap=my_cmap,
                     interpolation = 'none')
divider = make_axes_locatable(ax[0, 2])
cax3 = divider.append_axes("bottom", size=colorbarsize, pad=pad)
fig.colorbar(im, cax = cax3,
             label = "Measurement Time [s]",
             orientation="horizontal")

ax[0, 0].set_xlim(70, 130)
ax[0, 0].set_ylim(70, 130)
ax[0, 1].set_xlim(70, 130)
ax[0, 1].set_ylim(70, 130)
ax[0, 2].set_xlim(70, 130)
ax[0, 2].set_ylim(70, 130)

ax[0, 0].set_ylabel("[meters]")
ax[0, 0].set_xlabel("[meters]")
ax[0, 1].set_xlabel("[meters]")
ax[0, 2].set_xlabel("[meters]")

ax[0, 0].annotate("Source: Weapon", xy=(0.5, 1), xytext=(0, titlepad),
                  xycoords='axes fraction', textcoords='offset points',
                  size='large', ha='center', va='baseline')
ax[0, 1].annotate("Source: Cosmic", xy=(0.5, 1), xytext=(0, titlepad),
                  xycoords='axes fraction', textcoords='offset points',
                  size='large', ha='center', va='baseline')
ax[0, 2].annotate("Measurement Time", xy=(0.5, 1), xytext=(0, titlepad),
                  xycoords='axes fraction', textcoords='offset points',
                  size='large', ha='center', va='baseline')

#plt.tight_layout()
plt.savefig(os.path.join("plots-for-presentation", "soviet-lower-m{}.pdf".format(m)), transparent=True)
plt.savefig(os.path.join("plots-for-presentation", "soviet-lower-m{}.png".format(m)), transparent=True)

plt.show()

################################################################################
# Combined plot Soviet Storage for Presentation
m = 3
sostplotz1 = 2
sostplotz2 = 7
pad = 0.5
titlepad = 5
colorbarsize = "8%"

tmin = 1
tmax = 24 * 3600

factor = 0.8
fig, ax = plt.subplots(nrows = 1, ncols = 3, squeeze = False,
                       figsize = (12 * factor, 5 * factor),
                       sharey = True
)

vmin = 0.0001
vmax = 10 ** (math.ceil(math.log(sostshieldedcosmiccurrentdata.max(), 10)))

data = sostshieldedcurrentdata[sostplotz2]
im = ax[0, 0].imshow(data,
                     norm=LogNorm(vmin = vmin, vmax = vmax),
                     cmap=my_cmap,
                     interpolation = 'none')
divider = make_axes_locatable(ax[0, 0])
cax1 = divider.append_axes("bottom", size=colorbarsize, pad=pad)
fig.colorbar(im, cax = cax1,
             label = "Current into 1 m$^3$ [neutron/s]",
             orientation="horizontal")

data = sostshieldedcosmiccurrentdata[sostplotz2]
im = ax[0, 1].imshow(data,
                     norm=LogNorm(vmin = vmin, vmax = vmax),
                     cmap=my_cmap,
                     interpolation = 'none')
divider = make_axes_locatable(ax[0, 1])
cax2 = divider.append_axes("bottom", size=colorbarsize, pad=pad)
fig.colorbar(im, cax = cax2,
             label = "Current into 1 m$^3$ [neutron/s]",
             orientation="horizontal")

data = m ** 2 * sostshieldedcosmiccurrentdata[sostplotz2] / (sostshieldedcurrentdata[sostplotz2] ** 2)
data[data == np.inf] = 0
im = ax[0, 2].imshow(data,
                     norm=LogNorm(vmin = tmin, vmax = tmax),
                     cmap=my_cmap,
                     interpolation = 'none')
divider = make_axes_locatable(ax[0, 2])
cax3 = divider.append_axes("bottom", size=colorbarsize, pad=pad)
fig.colorbar(im, cax = cax3,
             label = "Measurement Time [s]",
             orientation="horizontal")

ax[0, 0].set_xlim(70, 130)
ax[0, 0].set_ylim(70, 130)
ax[0, 1].set_xlim(70, 130)
ax[0, 1].set_ylim(70, 130)
ax[0, 2].set_xlim(70, 130)
ax[0, 2].set_ylim(70, 130)

ax[0, 0].set_ylabel("[meters]")
ax[0, 0].set_xlabel("[meters]")
ax[0, 1].set_xlabel("[meters]")
ax[0, 2].set_xlabel("[meters]")

ax[0, 0].annotate("Source: Weapon", xy=(0.5, 1), xytext=(0, titlepad),
                  xycoords='axes fraction', textcoords='offset points',
                  size='large', ha='center', va='baseline')
ax[0, 1].annotate("Source: Cosmic", xy=(0.5, 1), xytext=(0, titlepad),
                  xycoords='axes fraction', textcoords='offset points',
                  size='large', ha='center', va='baseline')
ax[0, 2].annotate("Measurement Time", xy=(0.5, 1), xytext=(0, titlepad),
                  xycoords='axes fraction', textcoords='offset points',
                  size='large', ha='center', va='baseline')

#plt.tight_layout()
plt.savefig(os.path.join("plots-for-presentation", "soviet-upper-m{}.pdf".format(m)), transparent=True)
plt.savefig(os.path.join("plots-for-presentation", "soviet-upper-m{}.png".format(m)), transparent=True)

plt.show()
