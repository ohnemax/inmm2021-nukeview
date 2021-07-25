###############################################################################
# Some standard modules
import os
import sys
import copy
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
###############################################################################

import subprocess

scriptdir = os.path.dirname(os.path.realpath(__file__))
crytestcommand = os.path.join(scriptdir, "../cry-with-openmc/build/CRY-1.7-prefix/src/CRY-1.7/test/testMain")

defaultfile = """returnNeutrons 1
returnProtons 0
returnGammas 0
returnElectrons 0
returnMuons 0
returnPions 0
latitude 90.0
altitude 0
subboxLength 100
"""
width = 100

setupfile = "cry-check-distribution.file"
datadict = {
    "date": [],
    "npsm2": [],
    "npsm2-below-20MeV": [],
    "npsm2-below-200MeV": []
}
for year in range(2008, 2022):
    for month in range(1, 13, 3):
        print(year, month)
        with open(setupfile, 'w') as f:
            f.write(defaultfile)
            f.write("date {:d}-1-{:d}\n".format(month, year))
        f.close()
        
        proc = subprocess.Popen([crytestcommand, 'cry-check-distribution.file', '100000'],
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, erroroutput = proc.communicate()
        output = output.decode("utf-8") 
        lines = output.splitlines()

        xvalues = []
        yvalues = []
        evalues = []
        count = 0
        count20 = 0
        count200 = 0
        time = float(lines[-1].split(" ")[3])
        for line in lines:
            if line.startswith("Secondary"):
                e = float(line.split(" ")[3][3:])
                count += 1
                if e <= 200:
                    count200 += 1
                    if e <= 20:
                        count20 += 1
                evalues.append(e)
                xvalues.append(float(line.split(" ")[5]))
                yvalues.append(float(line.split(" ")[6]))
        datadict['date'].append("{:d}-{:02d}-01".format(year, month))
        datadict['npsm2'].append(count / time / (width ** 2))
        datadict['npsm2-below-20MeV'].append(count20 / time / (width ** 2))
        datadict['npsm2-below-200MeV'].append(count200 / time / (width ** 2))

df = pd.DataFrame(datadict)
df.to_csv("cry-check-distribution.csv")

# for last date
plt.plot(xvalues, yvalues, ".")
plt.title("source distribution")
plt.show()

plt.hist(xvalues, bins=100)
plt.title("distribution of x coordinates")
plt.show()

plt.hist(yvalues, bins=100)
plt.title("distribution of y coordinates")
plt.show()

plt.hist(evalues, bins=1000, range = (0, 100), log=True)
plt.xlabel("Energy / MeV")
plt.title("distribution of energy")
plt.show()

plt.hist(evalues, bins=1000, range = (0, 10), log=True)
plt.xlabel("Energy / MeV")
plt.title("distribution of energy")
plt.show()


