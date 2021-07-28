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
import argparse

if not os.path.exists("calculation.json"):
    print("Please execute script from inside of a simulation directory (contains calculation.json)")
    exit(-1)

with open("calculation.json", 'r') as f:
    calcsettings = json.load(f)
    ages = calcsettings['ages']
    geometries = calcsettings['geometries']
    particles = calcsettings['particles']
    f.close()

for geo in geometries:
    print(os.path.join("eigenvalue", geo))
    openmc.run(cwd = os.path.join("eigenvalue", geo))
    openmc.run(cwd = os.path.join("fixed-source", geo))
for age in ages:
    openmc.run(cwd = os.path.join("eigenvalue", "full-{:.2f}".format(age)))
    openmc.run(cwd = os.path.join("fixed-source", "full-{:.2f}".format(age)))    
