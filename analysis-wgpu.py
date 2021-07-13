###############################################################################
# Some standard modules
import os
import sys
import copy
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
###############################################################################

import openmc
import json

basepath = '20210609-100k'

with open(os.path.join(basepath, "calculation.json"), 'r') as f:
    calcsettings = json.load(f)
    ages = calcsettings['ages']
    geometries = calcsettings['geometries']
    f.close()
evpath = os.path.join(basepath, "eigenvalue")
fspath = os.path.join(basepath, "fixed-source")

radii = [5-0.75, 5,
         5 + 2,
         5 + 2 + 3,
         5 + 2 + 3 + 10,
         5 + 2 + 3 + 10 + 1,
         5 + 2 + 3 + 10 + 1 + 20]
volumes = [4 / 3 * np.pi * radii[0] ** 3]
for i in range(1, len(radii)):
    volumes.append(4 / 3 * np.pi * radii[i] ** 3)

batches = 120
agedfdict = {'type': [],
             'ages': [],
             'keff': [],
             'M': [],
             'kefftally': [],
             'Mtally': []}
fluxdf = pd.DataFrame(columns = ['type', 'cell', 'energy low [eV]', 'energy high [eV]', 'nuclide', 'score', 'mean', 'std. dev.', 'age'])
sourcedf = pd.DataFrame(columns = ['type', 'energy low [eV]', 'energy high [eV]', 'p', 'age'])

for age in ages:
    agedfdict['ages'].append(age)
    agedfdict['type'].append('eigenvalue')
    summary = openmc.Summary(os.path.join(evpath, "full-{:02d}/summary.h5".format(age)))
    sp = openmc.StatePoint(os.path.join(evpath, "full-{:02d}/statepoint.{:d}.h5".format(age, batches)))
    agedfdict['keff'].append(sp.k_combined)
    agedfdict['M'].append(1 / (1 - sp.k_combined))

    fiss_rate = sp.get_tally(name='fiss. rate')
    abs_rate = sp.get_tally(name='abs. rate')

    # Get the leakage tally
    leak = sp.get_tally(name='leakage')
    leak = leak.summation(filter_type=openmc.MeshSurfaceFilter, remove_filter=True)

    # Compute k-infinity using tally arithmetic
    keff = fiss_rate / (abs_rate + leak)
    keffval = keff.get_pandas_dataframe().loc[0, 'mean']
    agedfdict['kefftally'].append(keffval)
    agedfdict['Mtally'].append(1/(1-keffval))
    
    fluxtally = sp.get_tally(name='flux')
    tempdf = fluxtally.get_pandas_dataframe()
    tempdf['age'] = age
    tempdf['type'] = 'eigenvalue'
    fluxdf = pd.concat([fluxdf, tempdf])

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

    agedfdict['ages'].append(age)
    agedfdict['type'].append('fixed-source')
    summary = openmc.Summary(os.path.join(fspath, "full-{:02d}/summary.h5".format(age)))
    sp = openmc.StatePoint(os.path.join(fspath, "full-{:02d}/statepoint.{:d}.h5".format(age, batches)))
    agedfdict['keff'].append(np.nan)
    agedfdict['M'].append(np.nan)

    fiss_rate = sp.get_tally(name='fiss. rate')
    abs_rate = sp.get_tally(name='abs. rate')

    # Get the leakage tally
    leak = sp.get_tally(name='leakage')
    leak = leak.summation(filter_type=openmc.MeshSurfaceFilter, remove_filter=True)

    # Compute k-infinity using tally arithmetic
    keff = fiss_rate / (abs_rate + leak)
    keffval = keff.get_pandas_dataframe().loc[0, 'mean']
    agedfdict['kefftally'].append(keffval)
    agedfdict['Mtally'].append(1/(1-keffval))
    
    fluxtally = sp.get_tally(name='flux')
    tempdf = fluxtally.get_pandas_dataframe()
    tempdf['age'] = age
    tempdf['type'] = 'fixed-source'
    fluxdf = pd.concat([fluxdf, tempdf])

    tallyenergies = fluxtally.find_filter(openmc.EnergyFilter).values
    # sourcep, sourcebinedges = np.histogram(sp.source['E'], tallyenergies, density=True)
    #print(sum(sourcep * np.diff(tallyenergies)))
    sourcedfdict = {
        'energy low [eV]': [],
        'energy high [eV]': [],
        'p': []
    }
    tempdf = pd.DataFrame(sourcedfdict)
    tempdf['age'] = age
    sourcedf = pd.concat([sourcedf, tempdf])

    
agedf = pd.DataFrame(agedfdict)

print(fluxdf[(fluxdf['type'] == 'fixed-source') & (fluxdf['cell'] == 7)])

fluxdf['lethargy low'] = fluxdf['mean'] * fluxdf['energy low [eV]']
for i in range(1, 8):
    print(sum(fluxdf[fluxdf['cell'] == i]['mean']))
    fluxdf.loc[fluxdf['cell'] == i, 'mean'] = fluxdf.loc[fluxdf['cell'] == i, 'mean'] / volumes[i - 1]
    print(sum(fluxdf[fluxdf['cell'] == i]['mean']))

print(fluxdf[(fluxdf['type'] == 'fixed-source') & (fluxdf['cell'] == 7)])

##
fig, ax = plt.subplots()
fluxdfsel = fluxdf[(fluxdf['age'] == 0) & (fluxdf['type'] == 'fixed-source')]
for label, df in fluxdfsel.groupby('cell'):
    df.plot(x = 'energy low [eV]',
            y = 'mean',
            ax = ax,
            logx = True,
            logy = True,
            label = label)
plt.grid()
plt.show()
##
fig, ax = plt.subplots()
fluxdfsel = fluxdf[(fluxdf['age'] == 0) & (fluxdf['type'] == 'fixed-source')]
for label, df in fluxdfsel.groupby('cell'):
    df.plot(x = 'energy low [eV]',
            y = 'lethargy low',
            ax = ax,
            logx = True,
            logy = True,
            label = label)
plt.grid()
plt.savefig("lethargy-fixed-source.png")
plt.show()
##
fig, ax = plt.subplots()
for label, df in sourcedf.groupby('age'):
    df.plot(x = 'energy low [eV]',
            y = 'p',
            ax = ax,
            logx = True,
            label = label)
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

fluxdfsel = fluxdf[fluxdf['cell'] == 2]
for label, df in fluxdfsel.groupby('age'):
    df.plot(x = 'energy low [eV]',
            y = 'mean',
            ax = ax,
            logx = True,
            label = label)
plt.show()
