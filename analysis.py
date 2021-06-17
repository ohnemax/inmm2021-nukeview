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

basepath = '20210605-100k-evonly'

with open(os.path.join(basepath, "calculation.json"), 'r') as f:
    calcsettings = json.load(f)
    ages = calcsettings['ages']
    geometries = calcsettings['geometries']
    f.close()
evpath = os.path.join(basepath, "eigenvalue")
fspath = os.path.join(basepath, "fixed-source")

batches = 120
agedfdict = {'ages': [],
             'keff': [],
             'M': [],
             'kefftally': [],
             'Mtally': []}
fluxdf = pd.DataFrame(columns = ['cell', 'energy low [eV]', 'energy high [eV]', 'nuclide', 'score', 'mean', 'std. dev.', 'age'])
sourcedf = pd.DataFrame(columns = ['energy low [eV]', 'energy high [eV]', 'p', 'age'])

for age in ages:
    agedfdict['ages'].append(age)
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

    
agedf = pd.DataFrame(agedfdict)

fluxdf['lethargy low'] = fluxdf['mean'] * fluxdf['energy low [eV]']
##
fig, ax = plt.subplots()
fluxdfsel = fluxdf[fluxdf['age'] == 0]
for label, df in fluxdfsel.groupby('cell'):
    df.plot(x = 'energy low [eV]',
            y = 'mean',
            ax = ax,
            logx = True,
            label = label)
plt.show()
##
fig, ax = plt.subplots()
fluxdfsel = fluxdf[fluxdf['age'] == 0]
for label, df in fluxdfsel.groupby('cell'):
    df.plot(x = 'energy low [eV]',
            y = 'lethargy low',
            ax = ax,
            logx = True,
            logy = True,
            label = label)
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
