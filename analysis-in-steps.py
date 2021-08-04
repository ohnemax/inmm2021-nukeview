import openmc
import json
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib
import math
import copy
import sys

def getrectmeshtally(basepath, batches, name):
    sp = openmc.StatePoint(os.path.join(basepath, "statepoint.{:d}.h5".format(batches)))
    tally = sp.get_tally(name=name)
    xwidth = len(tally.filters[0].mesh.x_grid) - 1
    ywidth = len(tally.filters[0].mesh.y_grid) - 1
    zwidth = len(tally.filters[0].mesh.z_grid) - 1
    fluxdata = tally.get_reshaped_data()[:,0,0]
    fluxdata = fluxdata.reshape((zwidth, ywidth, xwidth))
    relerrordata = tally.get_reshaped_data(value="rel_err")[:,0,0]
    relerrordata = relerrordata.reshape((zwidth, ywidth, xwidth))
    np.nan_to_num(relerrordata, copy = False)
    return fluxdata, relerrordata

def getregmeshtally(basepath, batches, name):
    sp = openmc.StatePoint(os.path.join(basepath, "statepoint.{:d}.h5".format(batches)))
    tally = sp.get_tally(name=name)
    (xwidth, ywidth, zwidth) = tally.filters[0].mesh.dimension
    fluxdata = tally.get_reshaped_data()[:,0,0]
    fluxdata = fluxdata.reshape((zwidth, ywidth, xwidth))
    relerrordata = tally.get_reshaped_data(value="rel_err")[:,0,0]
    relerrordata = relerrordata.reshape((zwidth, ywidth, xwidth))
    np.nan_to_num(relerrordata, copy = False)
    return fluxdata, relerrordata

################################################################################
print("DIFFERENCE BETWEEN SURVIVAL BIASING ON & OFF")
print("Büchel, 100 million particles")


plotpath = "plots/survival-on-off-buechel-100M/"
batches = 10

basepath = "cluster-results/bust-20210712-run-100M"
fluxdata, relerrordata = getrectmeshtally(basepath, batches, "flux")

basepath = "cluster-results/bust-20210712-run-100M-survial-biasing"
survivalfluxdata, survivalrelerrordata = getrectmeshtally(basepath, batches, "flux")


# thanks to https://stackoverflow.com/questions/9455044/problems-with-zeros-in-matplotlib-colors-lognorm
my_cmap = copy.copy(matplotlib.cm.get_cmap('viridis')) # copy the default cmap
my_cmap.set_bad(my_cmap(0))

print("Total")
print("Non-zero cells default: {:d} of {:d}".format(len(fluxdata[fluxdata > 0]), len(fluxdata.flatten())))
print("Non-zero cells survival biasing: {:d} of {:d}".format(len(survivalfluxdata[survivalfluxdata > 0]), len(survivalfluxdata.flatten())))
print("Average relative error default: {:f}".format(sum(relerrordata.flatten()) / len(fluxdata[fluxdata > 0])))
print("Average relative error survival biasing: {:f}".format(sum(survivalrelerrordata.flatten()) / len(survivalfluxdata[survivalfluxdata > 0])))

for z in range(len(fluxdata)):
    print("z = {:d}m".format(z))
    print("Non-zero cells default: {:d} of {:d}".format(len(fluxdata[z][fluxdata[z] > 0]), len(fluxdata[z].flatten())))
    print("Non-zero cells survival biasing: {:d} of {:d}".format(len(survivalfluxdata[z][survivalfluxdata[z] > 0]), len(survivalfluxdata[z].flatten())))
    print("Average relative error default: {:f}".format(sum(relerrordata[z].flatten()) / len(fluxdata[z][fluxdata[z] > 0])))
    print("Average relative error survival biasing: {:f}".format(sum(survivalrelerrordata[z].flatten()) / len(survivalfluxdata[z][survivalfluxdata[z] > 0])))
    if not os.path.exists(os.path.join(plotpath, "compare-z-{:02d}.png".format(z))):
        fig, ax = plt.subplots(2, 3, figsize=(30, 10))
        im = ax[0][0].imshow(fluxdata[z], norm=LogNorm(vmin=0.000001, vmax=10), cmap=my_cmap)
        ax[0][0].title.set_text("default, max: {:.3e}".format(fluxdata[z].max()))
        fig.colorbar(im, ax = ax[0][0])
        
        im = ax[0][1].imshow(survivalfluxdata[z], norm=LogNorm(vmin=0.000001, vmax=10), cmap=my_cmap)
        ax[0][1].title.set_text("survival biasing, max: {:.3e}".format(survivalfluxdata[z].max()))
        fig.colorbar(im, ax = ax[0][1])
        
        difference = np.abs(fluxdata[z] - survivalfluxdata[z])
        im = ax[0][2].imshow(difference, norm=LogNorm(vmin=0.0000001, vmax=0.01), cmap=my_cmap)
        ax[0][2].title.set_text("difference, max: {:.3e}".format(difference.max()))
        fig.colorbar(im, ax = ax[0][2])

        im = ax[1][0].imshow(relerrordata[z], norm=LogNorm(vmin=0.000001, vmax=1), cmap=my_cmap)
        ax[1][0].title.set_text("default, rel. error")
        fig.colorbar(im, ax = ax[1][0])

        im = ax[1][1].imshow(survivalrelerrordata[z], norm=LogNorm(vmin=0.000001, vmax=1), cmap=my_cmap)
        ax[1][1].title.set_text("survival biasing, rel. error")
        fig.colorbar(im, ax = ax[1][1])

        fig.delaxes(ax[1][2])
        plt.savefig(os.path.join(plotpath, "compare-z-{:02d}.png".format(z)))
        plt.savefig(os.path.join(plotpath, "compare-z-{:02d}.pdf".format(z)))
        plt.close()

################################################################################
print("DIFFERENCE BETWEEN SURVIVAL BIASING ON & OFF")
print("Büchel, 10 million particles")


plotpath = "plots/survival-on-off-buechel-10M/"
batches = 10

basepath = "cluster-results/bust-20210712-run-10M"
fluxdata, relerrordata = getrectmeshtally(basepath, batches, "flux")

basepath = "cluster-results/bust-20210712-run-10M-survial-biasing"
survivalfluxdata, survivalrelerrordata = getrectmeshtally(basepath, batches, "flux")


# thanks to https://stackoverflow.com/questions/9455044/problems-with-zeros-in-matplotlib-colors-lognorm
my_cmap = copy.copy(matplotlib.cm.get_cmap('viridis')) # copy the default cmap
my_cmap.set_bad(my_cmap(0))

print("Total")
print("Non-zero cells default: {:d} of {:d}".format(len(fluxdata[fluxdata > 0]), len(fluxdata.flatten())))
print("Non-zero cells survival biasing: {:d} of {:d}".format(len(survivalfluxdata[survivalfluxdata > 0]), len(survivalfluxdata.flatten())))
print("Average relative error default: {:f}".format(sum(relerrordata.flatten()) / len(fluxdata[fluxdata > 0])))
print("Average relative error survival biasing: {:f}".format(sum(survivalrelerrordata.flatten()) / len(survivalfluxdata[survivalfluxdata > 0])))

for z in range(len(fluxdata)):
    print("z = {:d}m".format(z))
    print("Non-zero cells default: {:d} of {:d}".format(len(fluxdata[z][fluxdata[z] > 0]), len(fluxdata[z].flatten())))
    print("Non-zero cells survival biasing: {:d} of {:d}".format(len(survivalfluxdata[z][survivalfluxdata[z] > 0]), len(survivalfluxdata[z].flatten())))
    print("Average relative error default: {:f}".format(sum(relerrordata[z].flatten()) / len(fluxdata[z][fluxdata[z] > 0])))
    print("Average relative error survival biasing: {:f}".format(sum(survivalrelerrordata[z].flatten()) / len(survivalfluxdata[z][survivalfluxdata[z] > 0])))
    if not os.path.exists(os.path.join(plotpath, "compare-z-{:02d}.png".format(z))):
        fig, ax = plt.subplots(2, 3, figsize=(30, 10))
        im = ax[0][0].imshow(fluxdata[z], norm=LogNorm(vmin=0.000001, vmax=10), cmap=my_cmap)
        ax[0][0].title.set_text("default, max: {:.3e}".format(fluxdata[z].max()))
        fig.colorbar(im, ax = ax[0][0])
        
        im = ax[0][1].imshow(survivalfluxdata[z], norm=LogNorm(vmin=0.000001, vmax=10), cmap=my_cmap)
        ax[0][1].title.set_text("survival biasing, max: {:.3e}".format(survivalfluxdata[z].max()))
        fig.colorbar(im, ax = ax[0][1])
        
        difference = np.abs(fluxdata[z] - survivalfluxdata[z])
        im = ax[0][2].imshow(difference, norm=LogNorm(vmin=0.0000001, vmax=0.01), cmap=my_cmap)
        ax[0][2].title.set_text("difference, max: {:.3e}".format(difference.max()))
        fig.colorbar(im, ax = ax[0][2])

        im = ax[1][0].imshow(relerrordata[z], norm=LogNorm(vmin=0.000001, vmax=1), cmap=my_cmap)
        ax[1][0].title.set_text("default, rel. error")
        fig.colorbar(im, ax = ax[1][0])

        im = ax[1][1].imshow(survivalrelerrordata[z], norm=LogNorm(vmin=0.000001, vmax=1), cmap=my_cmap)
        ax[1][1].title.set_text("survival biasing, rel. error")
        fig.colorbar(im, ax = ax[1][1])

        fig.delaxes(ax[1][2])
        plt.savefig(os.path.join(plotpath, "compare-z-{:02d}.png".format(z)))
        plt.savefig(os.path.join(plotpath, "compare-z-{:02d}.pdf".format(z)))
        plt.close()

################################################################################
print("DIFFERENCE BETWEEN FETTER & POINT SOURCE")
print("Büchel, 20 million particles")
print("Survival Biasing On")


plotpath = "plots/fetter-vs-point-buechel-20M/"
batches = 10

basepath = "cluster-results/bust-20210714-fetter-survival-run-20M"
fluxdata, relerrordata = getrectmeshtally(basepath, batches, "flux")

basepath = "cluster-results/bust-20210714-point-survival-run-20M"
pointfluxdata, pointrelerrordata = getrectmeshtally(basepath, batches, "flux")


# thanks to https://stackoverflow.com/questions/9455044/problems-with-zeros-in-matplotlib-colors-lognorm
my_cmap = copy.copy(matplotlib.cm.get_cmap('viridis')) # copy the default cmap
my_cmap.set_bad(my_cmap(0))

print("Total")
print("Non-zero cells fetter model source: {:d} of {:d}".format(len(fluxdata[fluxdata > 0]), len(fluxdata.flatten())))
print("Non-zero cells point source: {:d} of {:d}".format(len(pointfluxdata[pointfluxdata > 0]), len(pointfluxdata.flatten())))
print("Average relative error fetter model source: {:f}".format(sum(relerrordata.flatten()) / len(fluxdata[fluxdata > 0])))
print("Average relative error point source: {:f}".format(sum(pointrelerrordata.flatten()) / len(pointfluxdata[pointfluxdata > 0])))

for z in range(len(fluxdata)):
    print("z = {:d}m".format(z))
    print("Non-zero cells fetter model source: {:d} of {:d}".format(len(fluxdata[z][fluxdata[z] > 0]), len(fluxdata[z].flatten())))
    print("Non-zero cells point source: {:d} of {:d}".format(len(pointfluxdata[z][pointfluxdata[z] > 0]), len(pointfluxdata[z].flatten())))
    print("Average relative error fetter model source: {:f}".format(sum(relerrordata[z].flatten()) / len(fluxdata[z][fluxdata[z] > 0])))
    print("Average relative error point source: {:f}".format(sum(pointrelerrordata[z].flatten()) / len(pointfluxdata[z][pointfluxdata[z] > 0])))
    if not os.path.exists(os.path.join(plotpath, "compare-z-{:02d}.png".format(z))):
        fig, ax = plt.subplots(2, 3, figsize=(30, 10))
        im = ax[0][0].imshow(fluxdata[z], norm=LogNorm(vmin=0.000001, vmax=10), cmap=my_cmap)
        ax[0][0].title.set_text("fetter model source, max: {:.3e}".format(fluxdata[z].max()))
        fig.colorbar(im, ax = ax[0][0])
        
        im = ax[0][1].imshow(pointfluxdata[z], norm=LogNorm(vmin=0.000001, vmax=10), cmap=my_cmap)
        ax[0][1].title.set_text("point source, max: {:.3e}".format(pointfluxdata[z].max()))
        fig.colorbar(im, ax = ax[0][1])
        
        difference = np.abs(fluxdata[z] - pointfluxdata[z])
        im = ax[0][2].imshow(difference, norm=LogNorm(vmin=0.0000001, vmax=0.01), cmap=my_cmap)
        ax[0][2].title.set_text("|difference|, max: {:.3e}".format(difference.max()))
        fig.colorbar(im, ax = ax[0][2])

        im = ax[1][0].imshow(relerrordata[z], norm=LogNorm(vmin=0.000001, vmax=1), cmap=my_cmap)
        ax[1][0].title.set_text("fetter model source, rel. error")
        fig.colorbar(im, ax = ax[1][0])

        im = ax[1][1].imshow(pointrelerrordata[z], norm=LogNorm(vmin=0.000001, vmax=1), cmap=my_cmap)
        ax[1][1].title.set_text("point source, rel. error")
        fig.colorbar(im, ax = ax[1][1])

        fig.delaxes(ax[1][2])
        plt.savefig(os.path.join(plotpath, "compare-z-{:02d}.png".format(z)))
        plt.savefig(os.path.join(plotpath, "compare-z-{:02d}.pdf".format(z)))
        plt.close()
    
    
################################################################################
print("DIFFERENCE BETWEEN RECTILINEAR AND REGULAR MESH")
print("Büchel, 20 million particles")
print("With Fetter Model source and survival biasing on")

plotpath = "plots/rectilinear-vs-regular-buechel-20M/"
batches = 10

basepath = "cluster-results/bust-20210718-fetter-survial-run-20M"
fluxdata, relerrordata = getrectmeshtally(basepath, batches, "flux")

regularfluxdata, regularrelerrordata = getregmeshtally(basepath, batches, "flux, regular mesh")

# thanks to https://stackoverflow.com/questions/9455044/problems-with-zeros-in-matplotlib-colors-lognorm
my_cmap = copy.copy(matplotlib.cm.get_cmap('viridis')) # copy the default cmap
my_cmap.set_bad(my_cmap(0))

print("Total")
print("Non-zero cells rectilinear mesh: {:d} of {:d}".format(len(fluxdata[fluxdata > 0]), len(fluxdata.flatten())))
print("Non-zero cells regular mesh: {:d} of {:d}".format(len(regularfluxdata[regularfluxdata > 0]), len(regularfluxdata.flatten())))
print("Average relative error rectilinear mesh: {:f}".format(sum(relerrordata.flatten()) / len(fluxdata[fluxdata > 0])))
print("Average relative error regular mesh: {:f}".format(sum(regularrelerrordata.flatten()) / len(regularfluxdata[regularfluxdata > 0])))

for z in range(len(fluxdata)):
    print("z = {:d}m".format(z))
    print("Non-zero cells rectilinear mesh: {:d} of {:d}".format(len(fluxdata[z][fluxdata[z] > 0]), len(fluxdata[z].flatten())))
    print("Non-zero cells regular mesh: {:d} of {:d}".format(len(regularfluxdata[z][regularfluxdata[z] > 0]), len(regularfluxdata[z].flatten())))
    print("Average relative error rectilinear mesh: {:f}".format(sum(relerrordata[z].flatten()) / len(fluxdata[z][fluxdata[z] > 0])))
    print("Average relative error regular mesh: {:f}".format(sum(regularrelerrordata[z].flatten()) / len(regularfluxdata[z][regularfluxdata[z] > 0])))
    if not os.path.exists(os.path.join(plotpath, "compare-z-{:02d}.png".format(z))):
        fig, ax = plt.subplots(2, 3, figsize=(30, 10))
        im = ax[0][0].imshow(fluxdata[z], norm=LogNorm(vmin=0.000001, vmax=10), cmap=my_cmap)
        ax[0][0].title.set_text("rectilinear mesh, max: {:.3e}".format(fluxdata[z].max()))
        fig.colorbar(im, ax = ax[0][0])
        
        im = ax[0][1].imshow(regularfluxdata[z], norm=LogNorm(vmin=0.000001, vmax=10), cmap=my_cmap)
        ax[0][1].title.set_text("regular mesh, max: {:.3e}".format(regularfluxdata[z].max()))
        fig.colorbar(im, ax = ax[0][1])
        
        difference = np.abs(fluxdata[z] - regularfluxdata[z])
        im = ax[0][2].imshow(difference, norm=LogNorm(vmin=0.0000001, vmax=0.01), cmap=my_cmap)
        ax[0][2].title.set_text("difference, max: {:.3e}".format(difference.max()))
        fig.colorbar(im, ax = ax[0][2])

        im = ax[1][0].imshow(relerrordata[z], norm=LogNorm(vmin=0.000001, vmax=1), cmap=my_cmap)
        ax[1][0].title.set_text("rectilinear mesh, rel. error")
        fig.colorbar(im, ax = ax[1][0])

        im = ax[1][1].imshow(regularrelerrordata[z], norm=LogNorm(vmin=0.000001, vmax=1), cmap=my_cmap)
        ax[1][1].title.set_text("regular mesh, rel. error")
        fig.colorbar(im, ax = ax[1][1])

        fig.delaxes(ax[1][2])
        plt.savefig(os.path.join(plotpath, "compare-z-{:02d}.png".format(z)))
        plt.savefig(os.path.join(plotpath, "compare-z-{:02d}.pdf".format(z)))
        plt.close()

################################################################################
print("DIFFERENCE BETWEEN SURVIVAL BIASING WEIGHT CUTOFF (0.25 / 0.05)")
print("Büchel, 50 million particles")


plotpath = "plots/weight-cutoff-0.05-0.25-buechel-50M/"
batches = 10

basepath = "cluster-results/bust-20210725-fetter-survival-50M"
fluxdata, relerrordata = getregmeshtally(basepath, batches, "flux, regular mesh")

basepath = "cluster-results/bust-20210725-fetter-survival-0.05-50M"
survivalfluxdata, survivalrelerrordata = getregmeshtally(basepath, batches, "flux, regular mesh")


# thanks to https://stackoverflow.com/questions/9455044/problems-with-zeros-in-matplotlib-colors-lognorm
my_cmap = copy.copy(matplotlib.cm.get_cmap('viridis')) # copy the default cmap
my_cmap.set_bad(my_cmap(0))

print("Total")
print("Non-zero cells default: {:d} of {:d}".format(len(fluxdata[fluxdata > 0]), len(fluxdata.flatten())))
print("Non-zero cells survival biasing: {:d} of {:d}".format(len(survivalfluxdata[survivalfluxdata > 0]), len(survivalfluxdata.flatten())))
print("Average relative error default: {:f}".format(sum(relerrordata.flatten()) / len(fluxdata[fluxdata > 0])))
print("Average relative error survival biasing: {:f}".format(sum(survivalrelerrordata.flatten()) / len(survivalfluxdata[survivalfluxdata > 0])))

for z in range(len(fluxdata)):
    print("z = {:d}m".format(z))
    print("Non-zero cells default: {:d} of {:d}".format(len(fluxdata[z][fluxdata[z] > 0]), len(fluxdata[z].flatten())))
    print("Non-zero cells survival biasing: {:d} of {:d}".format(len(survivalfluxdata[z][survivalfluxdata[z] > 0]), len(survivalfluxdata[z].flatten())))
    print("Average relative error default: {:f}".format(sum(relerrordata[z].flatten()) / len(fluxdata[z][fluxdata[z] > 0])))
    print("Average relative error survival biasing: {:f}".format(sum(survivalrelerrordata[z].flatten()) / len(survivalfluxdata[z][survivalfluxdata[z] > 0])))
    if not os.path.exists(os.path.join(plotpath, "compare-z-{:02d}.png".format(z))):
        fig, ax = plt.subplots(2, 3, figsize=(30, 10))
        im = ax[0][0].imshow(fluxdata[z], norm=LogNorm(vmin=0.000001, vmax=10), cmap=my_cmap)
        ax[0][0].title.set_text("default, max: {:.3e}".format(fluxdata[z].max()))
        fig.colorbar(im, ax = ax[0][0])
        
        im = ax[0][1].imshow(survivalfluxdata[z], norm=LogNorm(vmin=0.000001, vmax=10), cmap=my_cmap)
        ax[0][1].title.set_text("survival biasing, max: {:.3e}".format(survivalfluxdata[z].max()))
        fig.colorbar(im, ax = ax[0][1])
        
        difference = np.abs(fluxdata[z] - survivalfluxdata[z])
        im = ax[0][2].imshow(difference, norm=LogNorm(vmin=0.0000001, vmax=0.01), cmap=my_cmap)
        ax[0][2].title.set_text("difference, max: {:.3e}".format(difference.max()))
        fig.colorbar(im, ax = ax[0][2])

        im = ax[1][0].imshow(relerrordata[z], norm=LogNorm(vmin=0.000001, vmax=1), cmap=my_cmap)
        ax[1][0].title.set_text("default, rel. error")
        fig.colorbar(im, ax = ax[1][0])

        im = ax[1][1].imshow(survivalrelerrordata[z], norm=LogNorm(vmin=0.000001, vmax=1), cmap=my_cmap)
        ax[1][1].title.set_text("survival biasing, rel. error")
        fig.colorbar(im, ax = ax[1][1])

        fig.delaxes(ax[1][2])
        plt.savefig(os.path.join(plotpath, "compare-z-{:02d}.png".format(z)))
        plt.savefig(os.path.join(plotpath, "compare-z-{:02d}.pdf".format(z)))
        plt.close()
        
