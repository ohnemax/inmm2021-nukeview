################################################################################
# Combined plot
m = 3
plotz = 3

tmin = 1
tmax = 24 * 3600

vmin = 0.0001
vmax = 10 ** (math.ceil(math.log(shieldedcosmiccurrentdata.max(), 10)))

factor = 1
fig, ax = plt.subplots(nrows = 1, ncols = 3, squeeze = False,
                       sharey = True,
                       figsize = (12 * factor, 5 * factor)
)

data = shieldedcurrentdata[plotz]
im = ax[0, 0].imshow(data,
                     norm=LogNorm(vmin = vmin, vmax = vmax),
                     cmap=my_cmap,
                     interpolation = 'none')
ax[0, 0].set_xlim(70, 130)
ax[0, 0].set_ylim(70, 130)
#fig.colorbar(im, ax = [ax[0, 0]], location = 'top', label = "Inward Current [neutron/s]")
data = shieldedcosmiccurrentdata[plotz]
im = ax[0, 1].imshow(data,
                     norm=LogNorm(vmin = vmin, vmax = vmax),
                     cmap=my_cmap,
                     interpolation = 'none')
ax[0, 1].set_xlim(70, 130)
ax[0, 1].set_ylim(70, 130)
fig.colorbar(im, ax = ax[0, 1],
             orientation = "vertical",
             fraction=0.046, pad=0.04,
             label = "Inward Current [neutron/s]")

data = m ** 2 * shieldedcosmiccurrentdata[plotz] / (shieldedcurrentdata[plotz] ** 2)
data[data == np.inf] = 0
im = ax[0, 2].imshow(data,
                     norm=LogNorm(vmin = tmin, vmax = tmax),
                     cmap=my_cmap,
                     interpolation = 'none')
ax[0, 2].set_xlim(70, 130)
ax[0, 2].set_ylim(70, 130)
fig.colorbar(im, ax = ax[0, 2],
             orientation = "vertical",
             fraction=0.046, pad=0.04,
             label = "Measurement Time [s]")

ax[0, 0].set_ylabel("y-axis [m]")
ax[0, 0].set_xlabel("x-axis [m]")
ax[0, 1].set_xlabel("x-axis [m]")
ax[0, 2].set_xlabel("x-axis [m]")

fig.tight_layout()
plt.show()

################################################################################
# Combined plot
m = 3
plotz = 3

tmin = 1
tmax = 24 * 3600

vmin = 0.0001
vmax = 10 ** (math.ceil(math.log(shieldedcosmiccurrentdata.max(), 10)))

factor = 1
fig, ax = plt.subplots(nrows = 1, ncols = 3, squeeze = False,
                       sharey = True,
                       figsize = (12 * factor, 5 * factor))

data = shieldedcurrentdata[plotz]
im = ax[0, 0].imshow(data,
                     norm=LogNorm(vmin = vmin, vmax = vmax),
                     cmap=my_cmap,
                     interpolation = 'none')
ax[0, 0].set_xlim(70, 130)
ax[0, 0].set_ylim(70, 130)
fig.colorbar(im, ax = [ax[0, 0]], location = 'top', label = "Inward Current [neutron/s]")
data = shieldedcosmiccurrentdata[plotz]
im = ax[0, 1].imshow(data,
                     norm=LogNorm(vmin = vmin, vmax = vmax),
                     cmap=my_cmap,
                     interpolation = 'none')
ax[0, 1].set_xlim(70, 130)
ax[0, 1].set_ylim(70, 130)
fig.colorbar(im, ax = ax[0, 1], orientation = "horizontal", label = "Inward Current [neutron/s]")

data = m ** 2 * shieldedcosmiccurrentdata[plotz] / (shieldedcurrentdata[plotz] ** 2)
data[data == np.inf] = 0
im = ax[0, 2].imshow(data,
                     norm=LogNorm(vmin = tmin, vmax = tmax),
                     cmap=my_cmap,
                     interpolation = 'none')
ax[0, 2].set_xlim(70, 130)
ax[0, 2].set_ylim(70, 130)
fig.colorbar(im, ax = ax[0, 2], orientation = "horizontal", label = "Measurement Time [s]")

ax[0, 0].set_ylabel("y-axis [m]")
ax[0, 0].set_xlabel("x-axis [m]")
ax[0, 1].set_xlabel("x-axis [m]")
ax[0, 2].set_xlabel("x-axis [m]")

fig.tight_layout()
plt.show()

################################################################################
# Combined plot
m = 3
plotz = 3

tmin = 1
tmax = 24 * 3600

vmin = 0.0001
vmax = 10 ** (math.ceil(math.log(shieldedcosmiccurrentdata.max(), 10)))

factor = 1
fig, ax = plt.subplots(nrows = 2, ncols = 3, squeeze = False, figsize = (12 * factor, 5 * factor))

data = shieldedcurrentdata[plotz]
im = ax[0, 0].imshow(data,
                     norm=LogNorm(vmin = vmin, vmax = vmax),
                     cmap=my_cmap,
                     interpolation = 'none')
ax[0, 0].set_xlim(70, 130)
ax[0, 0].set_ylim(70, 130)
fig.colorbar(im, cax = ax[1, 0], orientation = "horizontal")
data = shieldedcosmiccurrentdata[plotz]
im = ax[0, 1].imshow(data,
                     norm=LogNorm(vmin = vmin, vmax = vmax),
                     cmap=my_cmap,
                     interpolation = 'none')
ax[0, 1].set_xlim(70, 130)
ax[0, 1].set_ylim(70, 130)
fig.colorbar(im, cax = ax[1, 1], orientation = "horizontal")

data = m ** 2 * shieldedcosmiccurrentdata[plotz] / (shieldedcurrentdata[plotz] ** 2)
data[data == np.inf] = 0
im = ax[0, 2].imshow(data,
                     norm=LogNorm(vmin = tmin, vmax = tmax),
                     cmap=my_cmap,
                     interpolation = 'none')
ax[0, 2].set_xlim(70, 130)
ax[0, 2].set_ylim(70, 130)
fig.colorbar(im, cax = ax[1, 2], orientation = "horizontal")

plt.show()
exit(-1)
