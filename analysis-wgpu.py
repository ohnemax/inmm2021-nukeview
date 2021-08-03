###############################################################################
# Some standard modules
import os
import sys
import copy
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
###############################################################################
import json
from uncertainties import ufloat

import openmc

import helper

basepath = 'cluster-results/fetter-20210728-100k'
if len(sys.argv) == 2:
    basepath = sys.argv[1]

with open(os.path.join(basepath, "calculation.json"), 'r') as f:
    calcsettings = json.load(f)
    ages = calcsettings['ages']
    geometries = calcsettings['geometries']
    batches = calcsettings['batches']
    f.close()
evpath = os.path.join(basepath, "eigenvalue")
fspath = os.path.join(basepath, "fixed-source")
pspath = os.path.join(basepath, "point-source")

radiallayerthickness = [5 - 0.75,
                        0.75,
                        2,
                        3,
                        10,
                        1,
                        20
                        ]
radii = []
last = 0
for rad in radiallayerthickness:
    last = last + rad
    radii.append(last)
volumes = [4 / 3 * np.pi * radii[0] ** 3]
for i in range(1, len(radii)):
    volumes.append(4 / 3 * np.pi * radii[i] ** 3)
print(volumes)

pumass = 4000
weaponage = 0

puvec = helper.createwgpuvec()
age = weaponage
puvec.createagedvector([age])
sourcesum = 0
for iso in puvec.pudf.index:
    sourcesum += puvec.puagedf.loc[(iso, age), 'sfneutrons-ufloat'] * puvec.puagedf.loc[(iso, age), 'wo']
nps = sourcesum * pumass
print("The WPu weapon source ({:.2f} years old, {:.1f}g) emits {:.1f} neutrons/s".format(age, pumass, nps))

################################################################################
# read point source results

if os.path.exists(os.path.join(pspath, "statepoint.{:d}.h5".format(batches))):
    summary = os.path.join(pspath, "summary.h5")
    sp = openmc.StatePoint(os.path.join(pspath, "statepoint.{:d}.h5".format(batches)))

    tally = sp.get_tally(name = 'energy')
    newsourcedf = tally.get_pandas_dataframe()

fluxdf = pd.DataFrame(columns = ['type', 'age', 'cell', 'energy low [eV]', 'energy high [eV]', 'nuclide', 'score', 'mean', 'std. dev.'])
sourcedf = pd.DataFrame(columns = ['type', 'energy low [eV]', 'energy high [eV]', 'p', 'age'])
currentdf = pd.DataFrame(columns = ['type', 'age', 'surface', 'energy low [eV]', 'energy high [eV]', 'nuclide', 'score', 'mean', 'std. dev.', 'cellfrom'])
plaincurrentdf = pd.DataFrame(columns = ['type', 'age', 'surface', 'cell', 'cellfrom', 'nuclide', 'score', 'mean', 'std. dev.'])

ratedf = pd.DataFrame()

print("Reading different geometries")
for geom in geometries:
    print(geom)
    if os.path.exists(os.path.join(evpath, "{:s}/statepoint.{:d}.h5".format(geom, batches + 20))):
        summary = os.path.join(evpath, "{:s}/summary.h5".format(geom))
        sp = openmc.StatePoint(os.path.join(evpath, "{:s}/statepoint.{:d}.h5".format(geom, batches + 20)))
 
        fiss_rate = sp.get_tally(name='fiss. rate')
        abs_rate = sp.get_tally(name='abs. rate')

        # Get the leakage tally
        leak = sp.get_tally(name='leakage')
        leak = leak.summation(filter_type=openmc.MeshSurfaceFilter, remove_filter=True)

        ratetally = sp.get_tally(name='various reaction rates')
        tempdf = ratetally.get_pandas_dataframe()
        tempdf = tempdf.pivot(index='nuclide', columns='score', values=['mean', 'std. dev.'])
        tempdf.reset_index(inplace = True)
        tempdf['age'] = 0
        tempdf['type'] = "eigenvalue"
        tempdf['geometry'] = geom
        tempdf['leak'] = ufloat(leak.mean.flatten()[0], leak.std_dev.flatten()[0])
        tempdf['k-combined'] = sp.k_combined
        ratedf = pd.concat([ratedf, tempdf])

        # Compute k-infinity using tally arithmetic
        keff = fiss_rate / (abs_rate + leak)
        keffval = keff.get_pandas_dataframe().loc[0, 'mean']

    if os.path.exists(os.path.join(fspath, "{:s}/statepoint.{:d}.h5".format(geom, batches))):
        summary = os.path.join(fspath, "{:s}/summary.h5".format(geom))
        sp = openmc.StatePoint(os.path.join(fspath, "{:s}/statepoint.{:d}.h5".format(geom, batches)))

        fiss_rate = sp.get_tally(name='fiss. rate')
        abs_rate = sp.get_tally(name='abs. rate')

        # Get the leakage tally
        leak = sp.get_tally(name='leakage')
        leak = leak.summation(filter_type=openmc.MeshSurfaceFilter, remove_filter=True)

        ratetally = sp.get_tally(name='various reaction rates')
        tempdf = ratetally.get_pandas_dataframe()
        tempdf = tempdf.pivot(index='nuclide', columns='score', values=['mean', 'std. dev.'])
        tempdf.reset_index(inplace = True)
        tempdf['age'] = 0
        tempdf['type'] = "fixed-source"
        tempdf['geometry'] = geom
        tempdf['leak'] = ufloat(leak.mean.flatten()[0], leak.std_dev.flatten()[0])
        tempdf['k-combined'] = np.nan
        ratedf = pd.concat([ratedf, tempdf])

        # Compute k-infinity using tally arithmetic
        keff = fiss_rate / (abs_rate + leak)
        keffval = keff.get_pandas_dataframe().loc[0, 'mean']
    
print("Reading different ages")
for age in ages:
    if os.path.exists(os.path.join(evpath, "full-{:.2f}/statepoint.{:d}.h5".format(age, batches + 20))):
        print(age, "eigenvalue")
        summary = openmc.Summary(os.path.join(evpath, "full-{:.2f}/summary.h5".format(age)))
        sp = openmc.StatePoint(os.path.join(evpath, "full-{:.2f}/statepoint.{:d}.h5".format(age, batches + 20)))

        fiss_rate = sp.get_tally(name='fiss. rate')
        abs_rate = sp.get_tally(name='abs. rate')

        # Get the leakage tally
        leak = sp.get_tally(name='leakage')
        leak = leak.summation(filter_type=openmc.MeshSurfaceFilter, remove_filter=True)

        ratetally = sp.get_tally(name='various reaction rates')
        tempdf = ratetally.get_pandas_dataframe()
        tempdf = tempdf.pivot(index='nuclide', columns='score', values=['mean', 'std. dev.'])
        tempdf.reset_index(inplace = True)
        tempdf['age'] = age
        tempdf['type'] = "eigenvalue"
        tempdf['geometry'] = "full (age)"
        tempdf['leak'] = ufloat(leak.mean.flatten()[0], leak.std_dev.flatten()[0])
        tempdf['k-combined'] = sp.k_combined
        ratedf = pd.concat([ratedf, tempdf])

        # Compute k-infinity using tally arithmetic
        keff = fiss_rate / (abs_rate + leak)
        keffval = keff.get_pandas_dataframe().loc[0, 'mean']

        fluxtally = sp.get_tally(name='flux')
        tempdf = fluxtally.get_pandas_dataframe()
        tempdf['age'] = age
        tempdf['type'] = 'eigenvalue'
        fluxdf = pd.concat([fluxdf, tempdf])

        currenttally = sp.get_tally(name="surface current")
        tempdf = currenttally.get_pandas_dataframe()
        tempdf['age'] = age
        tempdf['type'] = 'eigenvalue'
        currentdf = pd.concat([currentdf, tempdf])

        plaincurrenttally = sp.get_tally(name="plain surface current")
        tempdf = plaincurrenttally.get_pandas_dataframe()
        tempdf['age'] = age
        tempdf['type'] = 'eigenvalue'
        plaincurrentdf = pd.concat([plaincurrentdf, tempdf])

        plaincurrenttally = sp.get_tally(name="plain cell current")
        tempdf = plaincurrenttally.get_pandas_dataframe()
        tempdf['age'] = age
        tempdf['type'] = 'eigenvalue'
        plaincurrentdf = pd.concat([plaincurrentdf, tempdf])

        plaincurrenttally = sp.get_tally(name="plain from cell current")
        tempdf = plaincurrenttally.get_pandas_dataframe()
        tempdf['age'] = age
        tempdf['type'] = 'eigenvalue'
        plaincurrentdf = pd.concat([plaincurrentdf, tempdf])
        
        tallyenergies = fluxtally.find_filter(openmc.EnergyFilter).values
        sourcep, sourcebinedges = np.histogram(sp.source['E'], tallyenergies, density=True)
        #print(sum(sourcep * np.diff(tallyenergies)))
        sourcedfdict = {
            'energy low [eV]': sourcebinedges[:-1],
            'energy high [eV]': sourcebinedges[1:],
            'p': sourcep
        }
        tempdf = pd.DataFrame(sourcedfdict)
        tempdf['age'] = age
        sourcedf = pd.concat([sourcedf, tempdf])

    if os.path.exists(os.path.join(fspath, "full-{:.2f}/statepoint.{:d}.h5".format(age, batches))):
        print(age, "fixed-source")
        summary = openmc.Summary(os.path.join(fspath, "full-{:.2f}/summary.h5".format(age)))
        sp = openmc.StatePoint(os.path.join(fspath, "full-{:.2f}/statepoint.{:d}.h5".format(age, batches)))

        # Get the leakage tally
        leak = sp.get_tally(name='leakage')
        leak = leak.summation(filter_type=openmc.MeshSurfaceFilter, remove_filter=True)

        ratetally = sp.get_tally(name='various reaction rates')
        tempdf = ratetally.get_pandas_dataframe()
        tempdf = tempdf.pivot(index='nuclide', columns='score', values=['mean', 'std. dev.'])
        tempdf.reset_index(inplace = True)
        tempdf['age'] = age
        tempdf['type'] = "fixed-source"
        tempdf['geometry'] = "full (age)"
        tempdf['leak'] = ufloat(leak.mean.flatten()[0], leak.std_dev.flatten()[0])
        tempdf['k-combined'] = sp.k_combined
        ratedf = pd.concat([ratedf, tempdf])
        
        fiss_rate = sp.get_tally(name='fiss. rate')
        abs_rate = sp.get_tally(name='abs. rate')

        # Compute k-infinity using tally arithmetic
        keff = fiss_rate / (abs_rate + leak)
        keffval = keff.get_pandas_dataframe().loc[0, 'mean']

        fluxtally = sp.get_tally(name='flux')
        tempdf = fluxtally.get_pandas_dataframe()
        tempdf['age'] = age
        tempdf['type'] = 'fixed-source'
        fluxdf = pd.concat([fluxdf, tempdf])

        currenttally = sp.get_tally(name="surface current")
        tempdf = currenttally.get_pandas_dataframe()
        tempdf['age'] = age
        tempdf['type'] = 'fixed-source'
        currentdf = pd.concat([currentdf, tempdf])

        plaincurrenttally = sp.get_tally(name="plain surface current")
        tempdf = plaincurrenttally.get_pandas_dataframe()
        tempdf['age'] = age
        tempdf['type'] = 'fixed-source'
        plaincurrentdf = pd.concat([plaincurrentdf, tempdf])

        plaincurrenttally = sp.get_tally(name="plain cell current")
        tempdf = plaincurrenttally.get_pandas_dataframe()
        tempdf['age'] = age
        tempdf['type'] = 'fixed-source'
        plaincurrentdf = pd.concat([plaincurrentdf, tempdf])

        plaincurrenttally = sp.get_tally(name="plain from cell current")
        tempdf = plaincurrenttally.get_pandas_dataframe()
        tempdf['age'] = age
        tempdf['type'] = 'fixed-source'
        plaincurrentdf = pd.concat([plaincurrentdf, tempdf])

        # tallyenergies = fluxtally.find_filter(openmc.EnergyFilter).values
        # # sourcep, sourcebinedges = np.histogram(sp.source['E'], tallyenergies, density=True)
        # #print(sum(sourcep * np.diff(tallyenergies)))
        # sourcedfdict = {
        #     'energy low [eV]': [],
        #     'energy high [eV]': [],
        #     'p': []
        # }
        # tempdf = pd.DataFrame(sourcedfdict)
        # tempdf['age'] = age
        # sourcedf = pd.concat([sourcedf, tempdf])

scores = ['absorption', 'fission', 'scatter', 'total', 
          '(n,2nd)', '(n,2n)', '(n,3n)', '(n,na)', '(n,n3a)', '(n,2na)',
          '(n,3na)', '(n,np)', '(n,n2a)', '(n,2n2a)', '(n,nd)', '(n,nt)',
          '(n,n3He)', '(n,nd2a)', '(n,nt2a)', '(n,4n)', '(n,2np)',
          '(n,3np)', '(n,n2p)',
          '(n,gamma)', '(n,p)', '(n,t)', '(n,3He)', '(n,a)',
          '(n,2a)', '(n,3a)', '(n,2p)', '(n,pa)', '(n,t2a)', '(n,d2a)',
          '(n,pd)', '(n,pt)', '(n,da)',
          'nu-fission', 'nu-scatter'
]
for s in scores:
    mean = ratedf[('mean', s)].values
    std = ratedf[('std. dev.', s)].values
    ratedf[s] = [ufloat(x[0], x[1]) for x in zip(mean, std)]
    ratedf.drop(columns = [('mean', s), ('std. dev.', s)], inplace = True)

ratedf['nproduction'] = ratedf['nu-fission'] + ratedf['nu-scatter'] - ratedf['scatter']
ratedf['nproduction2'] = ratedf['nu-fission'] + ratedf['(n,2n)'] + 2 * ratedf['(n,3n)']

ratedf['keff'] = ratedf['nproduction'] / (ratedf['absorption'] + ratedf['leak'])
ratedf['fissionkeff'] = ratedf['nu-fission'] / (ratedf['absorption'] + ratedf['leak'])
ratedf['M'] = 1 / (1 - ratedf['keff'])
ratedf['fissionM'] = 1 / (1 - ratedf['fissionkeff'])
ratedf['p_L'] = ratedf['leak'] / ratedf['M']
ratedf['p_c'] = (ratedf['absorption'] - ratedf['fission']) / ratedf['M']
ratedf['p'] = ratedf['fission'] / ratedf['M']

ratedf['producingabsorption'] = ratedf['fission'] + ratedf['(n,2n)']
ratedf['totalnu'] = ratedf['nproduction'] / ratedf['producingabsorption']
ratedf['fissionnu'] = ratedf['nu-fission'] / ratedf['fission']
ratedf['capture'] = ratedf['absorption'] - ratedf['producingabsorption']

ratedf['tp_c'] = ratedf['capture'] / ratedf['M']
ratedf['tp'] = ratedf['producingabsorption'] / ratedf['M']

ratedf['serber-alpha'] = ratedf['capture'] / ratedf['fission']
ratedf['serber-fissionalpha'] = (ratedf['absorption'] - ratedf['fission']) / ratedf['fission']
ratedf['serber-M_L'] = 1 + ratedf['producingabsorption'] * (ratedf['totalnu'] - 1 - ratedf['serber-alpha'])
ratedf['serber-fissionM_L'] = 1 + ratedf['fission'] * (ratedf['fissionnu'] - 1 - ratedf['serber-fissionalpha'])
ratedf['reilly-M_L'] = (1 - ratedf['tp'] - ratedf['tp_c']) / (1 - ratedf['tp'] * ratedf['totalnu'])

################################################################################
ratedf.set_index(['type', 'geometry', 'age'], inplace = True)
################################################################################
print("********************************************************************************")
print("  keff = {}".format(ratedf.loc[('fixed-source', 'full', 0), 'keff'].values[0]))
print("     M = {}".format(ratedf.loc[('fixed-source', 'full', 0), 'M'].values[0]))
print("  leak = {}".format(ratedf.loc[('fixed-source', 'full', 0), 'leak'].values[0]))
print()
print("   source rate = {:.0f} neutrons / s".format(nps))
print("  leakage rate = {:.0f} neutrons / s".format(nps * ratedf.loc[('fixed-source', 'full', 0), 'leak'].values[0]))
print("   source rate = {:.2u} neutrons / s".format(nps))
print("  leakage rate = {:.2u} neutrons / s".format(nps * ratedf.loc[('fixed-source', 'full', 0), 'leak'].values[0]))
print("********************************************************************************")
################################################################################
# print(fluxdf[(fluxdf['type'] == 'fixed-source') & (fluxdf['cell'] == 7)])

fluxdf['mean per volume'] = 0
for i in range(1, 8):
    # print(sum(fluxdf[fluxdf['cell'] == i]['mean']))
    fluxdf.loc[fluxdf['cell'] == i, 'mean per volume'] = fluxdf.loc[fluxdf['cell'] == i, 'mean'] / volumes[i - 1]
    # print(sum(fluxdf[fluxdf['cell'] == i]['mean']))
fluxdf['lethargy low'] = fluxdf['mean'] * fluxdf['energy low [eV]']
fluxdf['lethargy low per volume'] = fluxdf['mean per volume'] * fluxdf['energy low [eV]']

# print(fluxdf[(fluxdf['type'] == 'fixed-source') & (fluxdf['cell'] == 7)])

fluxdf['cell name'] = ""
for cellid in fluxdf['cell'].unique():
    fluxdf.loc[fluxdf['cell'] == cellid, 'cell name'] = summary.geometry.get_all_cells()[cellid].name

################################################################################
# compare keff

print(ratedf[['type', 'geometry', 'age', 'k-combined', 'keff', 'fissionkeff', 'M', 'fissionM', 'leak']])

################################################################################
# Outward cell current for age = 0, and energy distribution in vacuum by age

factor = 0.7
fig, ax = plt.subplots(nrows = 1, ncols = 2,
                       figsize=(12 * factor, 4.8 * factor), squeeze = False)

currentdfsel = currentdf[(currentdf['age'] == 0) & (currentdf['type'] == 'fixed-source')]
tempdf = currentdfsel.groupby(['surface', 'cellfrom']).sum()
tempdf.reset_index(inplace = True)
surfaces = tempdf['surface'].unique()
surfaces = [summary.geometry.get_all_cells()[s].name for s in surfaces]
meanvalues = tempdf[tempdf['surface'] == tempdf['cellfrom']]['mean'].values
#drop center & vacuum
surfaces = surfaces[1:-1]
meanvalues = meanvalues[1:-1]

ax[0][0].bar(surfaces, meanvalues * nps.n)
# ax[0][0].title.set_text("Outward current through cell surface")
ax[0][0].axes.set_ylabel("current [neutrons / s]")
ax[0][0].axhline(nps.n, color = "orange", linewidth = 2)
# ax[0][0].axes.set_xlabel("from Cell")
ax[0][0].grid()
fig.autofmt_xdate()

fluxdfsel = fluxdf[(fluxdf['cell'] == 7) & (fluxdf['type'] == 'fixed-source') & (fluxdf['age'] == 0)]
xvalues = ((fluxdfsel['energy low [eV]'] + fluxdfsel['energy low [eV]']) / 2).values
plotvalues = fluxdfsel['mean'].values
plotvalues *= nps.n

sourcevalues = newsourcedf['mean'].values
factor = 10
sourcevalues = sourcevalues / sourcevalues.sum() * plotvalues.sum() / factor
# sourcevalues *= nps
# fluxdfsel.plot(x = 'energy low [eV]',
#                y = 'total mean per volume',
#                ax = ax[0][1],
#                logx = True,
#                logy = True,
#                legend = False)
ax[0][1].plot(xvalues, plotvalues)
ax[0][1].plot(xvalues, sourcevalues)
ax[0][1].set_xscale("log")
# ax[0][1].set_yscale("log")
ax[0][1].axes.set_xlabel("Energy [eV]")
ax[0][1].axes.set_ylabel("Flux [neutron * cm]")
# ax[0][1].title.set_text("Neutron flux (outside)")
ax[0][1].grid()    

fig.tight_layout()
plt.savefig(os.path.join("plots-for-paper", "outward-current-flux.png"))
plt.savefig(os.path.join("plots-for-paper", "outward-current-flux.pdf"))
plt.show()


################################################################################
# Understanding currents
currentdfsel = currentdf[(currentdf['age'] == 0) & (currentdf['type'] == 'fixed-source')]
tempdf = currentdfsel.groupby(['surface', 'cellfrom']).sum()
tempdf.reset_index(inplace = True)
tempdf['absmean'] = np.abs(tempdf['mean'])

print("#################### CellFromFilter ####################")
calccellfrom = tempdf.groupby(['cellfrom']).sum()
cellfrom = plaincurrentdf[~plaincurrentdf['cellfrom'].isna()]
cellfrom = cellfrom[(cellfrom['age'] == 0) & (cellfrom['type'] == 'fixed-source')]
print(calccellfrom)
print(cellfrom)

print("#################### SurfaceFilter ####################")
calcsurface = tempdf.groupby(['surface']).sum()
surface = plaincurrentdf[~plaincurrentdf['surface'].isna()]
surface = surface[(surface['age'] == 0) & (surface['type'] == 'fixed-source')]
print(calcsurface)
print(surface)

print("#################### CellFilter ####################")
calccelllist = []
for i in range(1, 8):
    intoout = tempdf[(tempdf['surface'] == i - 1) & (tempdf['cellfrom'] == i - 1)]['absmean']
    outtoin = tempdf[(tempdf['surface'] == i) & (tempdf['cellfrom'] == i + 1)]['absmean']
    both = pd.concat([intoout,outtoin])
    val = both.sum()
    calccelllist.append(val)
calccell = pd.DataFrame({"cell": range(1, 8), "absmean": calccelllist})
cell = plaincurrentdf[~plaincurrentdf['cell'].isna()]
cell = cell[(cell['age'] == 0) & (cell['type'] == 'fixed-source')]
print(calccell)
print(cell)

################################################################################
# Plot surface & cell current
fig, ax = plt.subplots(2, 2, figsize=(10, 10), sharey = True)

currentdfsel = currentdf[(currentdf['age'] == 0) & (currentdf['type'] == 'fixed-source')]
surfaces = currentdfsel['surface'].unique()
totaldf = currentdfsel.groupby(['surface']).sum()
meanvalues = totaldf['mean'].values

ax[0][0].bar(surfaces, meanvalues)
ax[0][0].title.set_text("Current through surfaces")

cellcurrentdf = plaincurrentdf[~plaincurrentdf['cell'].isna()]
cellcurrentdfsel = cellcurrentdf[(cellcurrentdf['age'] == 0) & (cellcurrentdf['type'] == 'fixed-source')]
cells = cellcurrentdfsel['cell'].unique()
totaldf = cellcurrentdfsel.groupby(['cell']).sum()
meanvalues = totaldf['mean'].values

ax[0][1].bar(cells, meanvalues)
ax[0][1].title.set_text("Current through cells")

currentdfsel = currentdf[(currentdf['age'] == 0) & (currentdf['type'] == 'fixed-source')]
tempdf = currentdfsel.groupby(['surface', 'cellfrom']).sum()
tempdf.reset_index(inplace = True)
surfaces = tempdf['surface'].unique()
meanvalues = tempdf[tempdf['surface'] == tempdf['cellfrom']]['mean'].values

ax[1][0].bar(surfaces, meanvalues)
ax[1][0].title.set_text("Outward current per surface")

currentdfsel = currentdf[(currentdf['age'] == 0) & (currentdf['type'] == 'fixed-source')]
tempdf = currentdfsel.groupby(['surface', 'cellfrom']).sum()
tempdf.reset_index(inplace = True)
surfaces = tempdf['surface'].unique()
meanvalues = -tempdf[tempdf['surface'] == tempdf['cellfrom'] - 1]['mean'].values
meanvalues = np.append(meanvalues, 0)
ax[1][1].bar(surfaces, meanvalues)
ax[1][1].title.set_text("Inward current per surface")


fig.autofmt_xdate()
ax[0][0].grid()
ax[0][1].grid()
ax[1][0].grid()
ax[1][1].grid()
plt.show()

##
fluxdfsel = fluxdf[(fluxdf['age'] == 0) & (fluxdf['type'] == 'fixed-source')]
fig, ax = plt.subplots(2, 2, figsize=(16, 16))
# for label, df in fluxdfsel.groupby('cell name'):
#     df.plot(x = 'energy low [eV]',
#             y = 'mean',
#             ax = ax[0][0],
#             logx = True,
#             logy = True,
#             label = label)
#     df.plot(x = 'energy low [eV]',
#             y = 'lethargy low',
#             ax = ax[1][0],
#             logx = True,
#             logy = True,
#             label = label)
for label, df in fluxdfsel.groupby('cell'):
    # if label <= 2:
    #     print("****", label, "****")
    #     print(df.iloc[0])
    df.plot(x = 'energy low [eV]',
            y = 'mean',
            ax = ax[0][0],
            logx = True,
            logy = True,
            label = summary.geometry.get_all_cells()[label].name)
    df.plot(x = 'energy low [eV]',
            y = 'lethargy low',
            ax = ax[1][0],
            logx = True,
            logy = True,
            label = summary.geometry.get_all_cells()[label].name)
    df.plot(x = 'energy low [eV]',
            y = 'mean per volume',
            ax = ax[0][1],
            logx = True,
            logy = True,
            label = summary.geometry.get_all_cells()[label].name)
    df.plot(x = 'energy low [eV]',
            y = 'lethargy low per volume',
            ax = ax[1][1],
            logx = True,
            logy = True,
            label = summary.geometry.get_all_cells()[label].name)
ax[0][0].title.set_text("OpenMC result")
ax[0][1].title.set_text("OpenMC result divided by cell volume")
ax[0][0].legend() # necessary to reorder legend!
ax[1][0].legend() # necessary to reorder legend!
ax[0][1].legend()
ax[1][1].legend()

ax[0][0].grid()
ax[1][0].grid()
ax[0][1].grid()
ax[1][1].grid()
plt.show()


fluxdfsel = fluxdf[(fluxdf['age'] == 0) & (fluxdf['type'] == 'fixed-source')]
cells = fluxdfsel['cell name'].unique()
totaldf = fluxdfsel.groupby(['cell']).sum()
meanvalues = totaldf['mean'].values
meanpervolumevalues = totaldf['mean per volume'].values

fig, ax = plt.subplots(1, 2, figsize=(10, 5))
ax[0].bar(cells, meanvalues)
ax[1].bar(cells, meanpervolumevalues)
degress = 45
fig.autofmt_xdate()
ax[0].grid()
ax[1].grid()
plt.show()

    
fig, ax = plt.subplots(2, 2, figsize=(16, 16))
for label, df in fluxdfsel.groupby('cell name'):
    df.plot(x = 'energy low [eV]',
            y = 'mean',
            ax = ax[0][0],
            logx = True,
            logy = True,
            label = label)
    df.plot(x = 'energy low [eV]',
            y = 'lethargy low',
            ax = ax[1][0],
            logx = True,
            logy = True,
            label = label)
    df.plot(x = 'energy low [eV]',
            y = 'mean per volume',
            ax = ax[0][1],
            logx = True,
            logy = True,
            label = label)
    df.plot(x = 'energy low [eV]',
            y = 'lethargy low per volume',
            ax = ax[1][1],
            logx = True,
            logy = True,
            label = label)
ax[0][0].legend()
ax[1][0].legend()
ax[0][1].legend()
ax[1][1].legend()
ax[0][0].grid()
ax[1][0].grid()
ax[0][1].grid()
ax[1][1].grid()
plt.show()

##
fig, ax = plt.subplots()
for label, df in sourcedf.groupby('age'):
    df.plot(x = 'energy low [eV]',
            y = 'p',
            ax = ax,
            logx = True,
            label = label)
ax.legend()
plt.title("Source particle energies")
plt.show()

##
fig, ax = plt.subplots()
for label, df in sourcedf.groupby('age'):
    df.plot(x = 'energy low [eV]',
            y = 'p',
            ax = ax,
            logx = True,
            label = label)
ax.legend()

# fluxdfsel = fluxdf[fluxdf['cell'] == 2]
# for label, df in fluxdfsel.groupby('age'):
#     df.plot(x = 'energy low [eV]',
#             y = 'mean',
#             ax = ax,
#             logx = True,
#             label = label)
# ax.legend()
plt.show()
