# INMM 2021 "Fetter Model Revisited" - Inputs & Analysis

This repository contains scripts to create input files and read results used for the article "Fetter Model Revisited: Detecting Nuclear Weapons 30 Years Later" by (@ohnemax)[https://github.com/ohnemax/], (@jelfes)[https://github.com/jelfes] and (@cfichtlscherer)[https://github.com/cfichtlscherer/]. The article was presented at the INMM \& ESARDA Joint Virtual Annual Meeting\\ August 23-26 \& August 30-September 1, 2021.

All simulations use [OpenMC](https://docs.openmc.org/en/stable/), an open source Monte Carlo particle transport simulation. Certain parts also use [CRY](https://nuclear.llnl.gov/simulation/) as a special particle source [cry-with-openmc](https://github.com/ohnemax/cry-with-openmc). To carry out the simulations, the `cry-with-openmc` repository should be present in the same parent folder as this repository, and the library should be compiled.

Different sets of input files can be created using `input-<...>.py` scripts. Most scripts can be controlled through command line arguments, `script.py -h` provides help. After creating input files, `openmc` should be executed in the input directories. No results are included here due to size concerns. 

Reults can be analyzed with various analysis scripts, again correcponding to the input file names. These read out tally data and produce plots. For the article above, `analysis-wgpu.py`, `analysis-simple-cosmic-ray-neutrons.py`, and `analysis-buechel-soviet.py` were used.



