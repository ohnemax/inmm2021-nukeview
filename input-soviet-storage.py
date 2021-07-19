import openmc
import json
import os
import matplotlib.pyplot as plt
import math
import numpy as np
import sys
import argparse

import helper

cspath = "/openmc/openmc-data/v0.12/lib80x_hdf5/cross_sections.xml"

basepath = "sost-20210719-point"
fettersource = False
particles = 1
batches = 10
plot = False
weaponage = 0
survival = False

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--fetter", action="store_true",
                    help="add fetter model as source")
parser.add_argument("-p", "--plot", action="store_true",
                    help="plot some cuts through the geometry")
parser.add_argument("-s", "--survival", action="store_true",
                    help="turn on survival biasing")
parser.add_argument("-b", "--basepath", 
                    help="set basepath for calculation")
args = parser.parse_args()
if args.fetter:
    fettersource = True
if args.plot:
    plot = True
if args.survival:
    survival = True
if not args.basepath is None:
    basepath = args.basepath

if not os.path.exists(basepath):
    os.mkdir(basepath)

# dump basic settings to json file
with open(os.path.join(basepath, "calculation.json"), 'w') as f:
    json.dump({'particles': particles,
               'basepath': basepath,
               'cspath': cspath,
               'batches': batches,
               'fettersource': fettersource,
               'weaponage': weaponage},
              f)
    f.close()

###############################################################################
# Materials
###############################################################################

steelMat = openmc.Material(name = "Steel (for doors)")
steelMat.set_density("g/cm3", 7.82) # PNNL Compendium page 281
steelMat.add_element("C", 0.022831, "ao")
steelMat.add_element("Fe", 0.977169, "ao")

concreteRegularMat = openmc.Material(name = "Regular Concrete")
concreteRegularMat.set_density("g/cm3", 2.3) # PNNL Compendium page 112
concreteRegularMat.add_element("H", 0.16038, "ao")
concreteRegularMat.add_element("O", 0.563183, "ao")
concreteRegularMat.add_element("Na", 0.021365, "ao")
concreteRegularMat.add_element("Al", 0.021343, "ao")
concreteRegularMat.add_element("Si", 0.203231, "ao")
concreteRegularMat.add_element("Ca", 0.018595, "ao")
concreteRegularMat.add_element("Fe", 0.004246, "ao")

# 300kg carbon steel / m3
# https://www.yourspreadsheets.co.uk/reinforcement-estimates.html
steelweightperm3 = 300
steelvolumefraction = 300 / (7.82 * 1e4)
concretevolumefraction = 1 - steelvolumefraction
wallMat = openmc.Material.mix_materials([concreteRegularMat, steelMat],
                                        [concretevolumefraction, steelvolumefraction], 'vo')
wallMat.name = "Wall Material"

airMat = openmc.Material(name='Air, Dry (near sea level)')
airMat.set_density('g/cm3', 1.205e-3) # PNNL Compendium page 18
airMat.add_element('C', 0.000124, 'wo')
airMat.add_element('N', 0.755268, 'wo')
airMat.add_element('O', 0.231781, 'wo')
airMat.add_element('Ar', 0.012827, 'wo')

soilMat = openmc.Material(name='Earth, Typical Western US')
soilMat.set_density('g/cm3', 1.52) # PNNL Compendium p.118 (references Brewer 2009, which has a density of 5.5g/cm3)
soilMat.add_element('H', 0.316855, 'ao')
soilMat.add_element('O', 0.501581, 'ao')
soilMat.add_element('Al', 0.039951, 'ao')
soilMat.add_element('Si', 0.141613, 'ao')

#*******************************************************************************
# fetter model

# need pit radius already here
pitOR = 5
pitThickness = 0.75
pitIR = pitOR - pitThickness

puvec = helper.puvector()
puvec.setwo('Pu238', 0.00005)
puvec.setwo('Pu239', 0.993)
puvec.setwo('Pu240', 0.060)
puvec.setwo('Pu241', 0.0044)
puvec.setwo('Pu242', 0.00015)
puvec.setwo('Am241', 0)
ages = [weaponage]
puvec.createagedvector(ages)

pudf = puvec.pudf
puagedf = puvec.puagedf

wpuRho = 4000 / (4/3 * np.pi * (pitOR ** 3 - pitIR ** 3)) # density based on Fetter Model volume / mass
wpuMat = {}
for age in ages:
    wpuMat[age] = openmc.Material(name = "Weapon-grade Plutonium - {:02d} years old".format(age))
    wpuMat[age].set_density("g/cm3", wpuRho)
    wpuMat[age].add_nuclide("Pu238", puagedf.loc[('Pu238', age), 'wo'], "wo")
    wpuMat[age].add_nuclide("Pu239", puagedf.loc[('Pu239', age), 'wo'], "wo")
    wpuMat[age].add_nuclide("Pu240", puagedf.loc[('Pu240', age), 'wo'], "wo") #should be wo all along
    wpuMat[age].add_nuclide("Pu241", puagedf.loc[('Pu241', age), 'wo'], "wo") # need to consider for decay!
    wpuMat[age].add_nuclide("Pu242", puagedf.loc[('Pu242', age), 'wo'], "wo")
    wpuMat[age].add_nuclide("Am241", puagedf.loc[('Am241', age), 'wo'], "wo")
    wpuMat[age].add_element("O", 0.002, "wo")

berMat = openmc.Material(name = "Beryllium Reflector")
berMat.set_density("g/cm3", 1.848) #PNNL Compendium page 37
berMat.add_element("Be", 1)

tunMat = openmc.Material(name = "Tungsten Tamper")
tunMat.set_density("g/cm3", 19.3) #PNNL Compendium page 318
tunMat.add_element("W", 1)

hmxMat = openmc.Material(name = "HMX explosive")
hmxMat.set_density("g/cm3", 1.890) #PNNL Compendium page 125
hmxMat.add_element("H", 0.285714, "ao")
hmxMat.add_element("C", 0.142857, "ao")
hmxMat.add_element("N", 0.285714, "ao")
hmxMat.add_element("O", 0.285714, "ao")

aluMat = openmc.Material(name = "Aluminum Case")
aluMat.set_density("g/cm3", 2.6989) #PNNL Compendium page 20
aluMat.add_element("Al", 1)

materiallist = [steelMat, concreteRegularMat, wallMat, airMat, soilMat]
if fettersource:
    materiallist = materiallist + [wpuMat[weaponage], berMat, tunMat, hmxMat, aluMat]

materials = openmc.Materials(materiallist)
materials.cross_sections = cspath
materials.export_to_xml(os.path.join(basepath, "materials.xml"))

###############################################################################
# Geometry
###############################################################################

horizontalpixel = 175 / 26.4 # from drawing - person with 1.75 is 26.4px tall
verticalpixel = 175 / 26.4   # from drawing - person with 1.75 is 26.4px tall
thickA = 22.7 * horizontalpixel
thickB = 9.5 * horizontalpixel
thickC = 3.4 * verticalpixel
thickD = 10 # Main doors directly in cm
thickE = 5
thickF = 20
dista = 0
distb = 107.6 * horizontalpixel
distc = distb / 3
distd = 0
diste = 454 * horizontalpixel
distf = 48.6 * verticalpixel
distg = 36.1 * verticalpixel
disth = 25 * verticalpixel #dirt height - an assumption
disti = 0
distj = 0
distk = 0
distl = 100 * horizontalpixel #TODO measure
distm = 150 * horizontalpixel #TODO measure
distn = 25 * horizontalpixel #TODO measure
disto = 25 * horizontalpixel #TODO measure
distp = 500 # Directly in cm

surfaces = {}

# x surfaces for lower level (main)
surfaceincrements = {1: 0,
                     5: distp,
                     7: disth,
                     10: thickA,
                     20: distl,
                     30: thickA,
                     40: distb / 3,
                     50: distb / 3,
                     60: distb / 3,
                     70: thickB,
                     80: distb / 3,
                     90: distb / 3,
                     100: distb / 3,
                     110: thickB,
                     120: distb / 3,
                     130: distb / 3,
                     140: distb / 3,
                     150: thickB,
                     160: distb / 3,
                     170: distb / 3,
                     180: distb / 3,
                     190: thickA,
                     200: distl,
                     210: thickA,
                     215: disth,
                     220: distp
}
xwidth = sum(surfaceincrements.values()) 

val = 0
for sid, inc in surfaceincrements.items():
    val += inc
    surfaces[sid] = openmc.XPlane(val, surface_id = sid)

# Define outer surfaces with vacuum boundary type
surfaces[1].boundary_type = 'vacuum'
surfaces[220].boundary_type = 'vacuum'

# Special surfaces for doors
surfaces[2] = openmc.XPlane(surfaces[5].x0 + thickD, surface_id = 2)
surfaces[28] = openmc.XPlane(surfaces[30].x0 - thickD, surface_id = 28)
surfaces[182] = openmc.XPlane(surfaces[180].x0 + thickD, surface_id = 182)
surfaces[208] = openmc.XPlane(surfaces[210].x0 - thickD, surface_id = 208)
    
surfaces[37] = openmc.XPlane(surfaces[40].x0 - thickF, surface_id = 37)    
surfaces[53] = openmc.XPlane(surfaces[50].x0 + thickF, surface_id = 53)    
surfaces[77] = openmc.XPlane(surfaces[80].x0 - thickF, surface_id = 77)    
surfaces[93] = openmc.XPlane(surfaces[90].x0 + thickF, surface_id = 93)    
surfaces[117] = openmc.XPlane(surfaces[120].x0 - thickF, surface_id = 117)    
surfaces[133] = openmc.XPlane(surfaces[130].x0 + thickF, surface_id = 133)    
surfaces[157] = openmc.XPlane(surfaces[160].x0 - thickF, surface_id = 157)    
surfaces[173] = openmc.XPlane(surfaces[170].x0 + thickF, surface_id = 173)    

# surfaces for outer soil
surfaces[15] = openmc.XPlane(surfaces[20].x0 - disth, surface_id = 15)
surfaces[195] = openmc.XPlane(surfaces[190].x0 + disth, surface_id = 195)

# x surfaces for lower level (auxiliary)

# y surfaces for lower level (main)
surfaceincrements = {990: 0,
                     995: distp,
                     1000: disth,
                     1010: thickA,
                     1020: diste,
                     1030: thickA,
                     1040: distn,
                     1050: disto,
                     1060: thickA,
                     1070: disto,
                     1080: distn,
                     1090: thickA,
                     1100: distm,
                     1110: thickA,
                     1115: disth,
                     1120: distp
}
ywidth = sum(surfaceincrements.values()) 

val = 0
for sid, inc in surfaceincrements.items():
    val += inc
    surfaces[sid] = openmc.YPlane(val, surface_id = sid)

# Define outer surfaces with vacuum boundary type
surfaces[990].boundary_type = 'vacuum'
surfaces[1120].boundary_type = 'vacuum'
    

# surface for weapon doors
surfaces[1032] = openmc.YPlane(surfaces[1030].y0 + thickE, surface_id = 1032)

# y surfaces for lower level (auxiliary)



# z surfaces (main)
surfaceincrements = {2000: 0,
                     2010: thickA,
                     2020: distf,
                     2030: thickC,
                     2040: thickA - thickC,
                     2050: disth - (thickA - thickC),
                     2060: distg - disth,
                     2070: thickA,
                     2080: disth,
                     2090: distp
}
zwidth = sum(surfaceincrements.values())

val = 0
for sid, inc in surfaceincrements.items():
    val += inc
    surfaces[sid] = openmc.ZPlane(val, surface_id = sid)

# Define outer surfaces with vacuum boundary type
surfaces[2000].boundary_type = 'vacuum'
surfaces[2090].boundary_type = 'vacuum'

a, b, c, d = helper.planeparameter3points([surfaces[20].x0, surfaces[1000].y0, surfaces[2040].z0 + disth],
                                   [surfaces[30].x0, surfaces[1000].y0, surfaces[2040].z0 + disth],
                                   [surfaces[20].x0, surfaces[1020].y0, surfaces[2070].z0 + disth])
surfaces[2200] = openmc.Plane(a, b, c, d, surface_id = 2200)

a, b, c, d = helper.planeparameter3points([surfaces[20].x0, surfaces[1090].y0, surfaces[2070].z0 + disth],
                                   [surfaces[30].x0, surfaces[1090].y0, surfaces[2070].z0 + disth],
                                   [surfaces[20].x0, surfaces[1110].y0, surfaces[2040].z0 + disth])
surfaces[2210] = openmc.Plane(a, b, c, d, surface_id = 2210)

a, b, c, d = helper.planeparameter3points([surfaces[20].x0, surfaces[1000].y0 - (surfaces[2040].z0 + disth), 0],
                                   [surfaces[30].x0, surfaces[1000].y0 - (surfaces[2040].z0 + disth), 0],
                                   [surfaces[20].x0, surfaces[1000].y0, surfaces[2040].z0 + disth])
surfaces[2220] = openmc.Plane(a, b, c, d, surface_id = 2220)

a, b, c, d = helper.planeparameter3points([surfaces[20].x0, surfaces[1110].y0 + (surfaces[2040].z0 + disth), 0],
                                   [surfaces[30].x0, surfaces[1110].y0 + (surfaces[2040].z0 + disth), 0],
                                   [surfaces[20].x0, surfaces[1110].y0, surfaces[2040].z0 + disth])
surfaces[2230] = openmc.Plane(a, b, c, d, surface_id = 2230)


# cells: base concrete
baseconcretecells = []
baseconcretecells.append(openmc.Cell())
baseconcretecells[-1].region = +surfaces[2000] & -surfaces[2010] & +surfaces[5] & -surfaces[20] & +surfaces[1020] & -surfaces[1090]
baseconcretecells[-1].name = "Base Concrete Wing I"

baseconcretecells.append(openmc.Cell())
baseconcretecells[-1].region = +surfaces[2000] & -surfaces[2010] & +surfaces[20] & -surfaces[190] & +surfaces[1000] & -surfaces[1110]
baseconcretecells[-1].name = "Base Concrete Main"

baseconcretecells.append(openmc.Cell())
baseconcretecells[-1].region = +surfaces[2000] & -surfaces[2010] & +surfaces[190] & -surfaces[210] & +surfaces[1020] & -surfaces[1090]
baseconcretecells[-1].name = "Base Concrete Wing II"

# cells: lower level concrete walls
lowerlevelwallcells = []
lowerlevelwallcells.append(openmc.Cell())
lowerlevelwallcells[-1].region = +surfaces[2010] & -surfaces[2020] & +surfaces[5] & -surfaces[10] & +surfaces[1020] & -surfaces[1090]
lowerlevelwallcells[-1].name = "Wing I Wall low"

lowerlevelwallcells.append(openmc.Cell())
lowerlevelwallcells[-1].region = +surfaces[2010] & -surfaces[2020] & +surfaces[10] & -surfaces[20] & +surfaces[1080] & -surfaces[1090]
lowerlevelwallcells[-1].name = "Wing I Wall auxiliary side"

lowerlevelwallcells.append(openmc.Cell())
lowerlevelwallcells[-1].region = +surfaces[2010] & -surfaces[2020] & +surfaces[10] & -surfaces[20] & +surfaces[1020] & -surfaces[1030]
lowerlevelwallcells[-1].name = "Wing I Wall weapon side"

lowerlevelwallcells.append(openmc.Cell())
lowerlevelwallcells[-1].region = +surfaces[2010] & -surfaces[2020] & +surfaces[10] & -surfaces[30] & +surfaces[1050] & -surfaces[1060]
lowerlevelwallcells[-1].name = "Wing I Wall center"

lowerlevelwallcells.append(openmc.Cell())
lowerlevelwallcells[-1].region = +surfaces[2010] & -surfaces[2020] & +surfaces[20] & -surfaces[30] & +surfaces[1000] & -surfaces[1040]
lowerlevelwallcells[-1].name = "Weapon Wall low"

lowerlevelwallcells.append(openmc.Cell())
lowerlevelwallcells[-1].region = +surfaces[2010] & -surfaces[2020] & +surfaces[30] & -surfaces[180] & +surfaces[1000] & -surfaces[1010]
lowerlevelwallcells[-1].name = "Weapon Wall outer"

lowerlevelwallcells.append(openmc.Cell())
lowerlevelwallcells[-1].region = +surfaces[2010] & -surfaces[2020] & +surfaces[30] & -surfaces[40] & +surfaces[1020] & -surfaces[1030]
lowerlevelwallcells[-1].name = "Weapon Wall inner I"

lowerlevelwallcells.append(openmc.Cell())
lowerlevelwallcells[-1].region = +surfaces[2010] & -surfaces[2020] & +surfaces[50] & -surfaces[80] & +surfaces[1020] & -surfaces[1030]
lowerlevelwallcells[-1].name = "Weapon Wall inner II"

lowerlevelwallcells.append(openmc.Cell())
lowerlevelwallcells[-1].region = +surfaces[2010] & -surfaces[2020] & +surfaces[90] & -surfaces[120] & +surfaces[1020] & -surfaces[1030]
lowerlevelwallcells[-1].name = "Weapon Wall inner III"

lowerlevelwallcells.append(openmc.Cell())
lowerlevelwallcells[-1].region = +surfaces[2010] & -surfaces[2020] & +surfaces[130] & -surfaces[160] & +surfaces[1020] & -surfaces[1030]
lowerlevelwallcells[-1].name = "Weapon Wall inner IV"

lowerlevelwallcells.append(openmc.Cell())
lowerlevelwallcells[-1].region = +surfaces[2010] & -surfaces[2020] & +surfaces[170] & -surfaces[180] & +surfaces[1020] & -surfaces[1030]
lowerlevelwallcells[-1].name = "Weapon Wall inner V"

lowerlevelwallcells.append(openmc.Cell())
lowerlevelwallcells[-1].region = +surfaces[2010] & -surfaces[2020] & +surfaces[180] & -surfaces[190] & +surfaces[1000] & -surfaces[1040]
lowerlevelwallcells[-1].name = "Weapon Wall high"

lowerlevelwallcells.append(openmc.Cell())
lowerlevelwallcells[-1].region = +surfaces[2010] & -surfaces[2020] & +surfaces[60] & -surfaces[70] & +surfaces[1010] & -surfaces[1020]
lowerlevelwallcells[-1].name = "Weapon Wall thin low"

lowerlevelwallcells.append(openmc.Cell())
lowerlevelwallcells[-1].region = +surfaces[2010] & -surfaces[2020] & +surfaces[100] & -surfaces[110] & +surfaces[1010] & -surfaces[1020]
lowerlevelwallcells[-1].name = "Weapon Wall thin middle"

lowerlevelwallcells.append(openmc.Cell())
lowerlevelwallcells[-1].region = +surfaces[2010] & -surfaces[2020] & +surfaces[140] & -surfaces[150] & +surfaces[1010] & -surfaces[1020]
lowerlevelwallcells[-1].name = "Weapon Wall thin high"

lowerlevelwallcells.append(openmc.Cell())
lowerlevelwallcells[-1].region = +surfaces[2010] & -surfaces[2020] & +surfaces[20] & -surfaces[30] & +surfaces[1070] & -surfaces[1110]
lowerlevelwallcells[-1].name = "Auxiliary Wall low"

lowerlevelwallcells.append(openmc.Cell())
lowerlevelwallcells[-1].region = +surfaces[2010] & -surfaces[2020] & +surfaces[30] & -surfaces[180] & +surfaces[1100] & -surfaces[1110]
lowerlevelwallcells[-1].name = "Auxiliary Wall outer"

# add auxiliary doors
lowerlevelwallcells.append(openmc.Cell())
lowerlevelwallcells[-1].region = +surfaces[2010] & -surfaces[2020] & +surfaces[30] & -surfaces[180] & +surfaces[1080] & -surfaces[1090]
lowerlevelwallcells[-1].name = "Auxiliary Wall inner"

lowerlevelwallcells.append(openmc.Cell())
lowerlevelwallcells[-1].region = +surfaces[2010] & -surfaces[2020] & +surfaces[180] & -surfaces[190] & +surfaces[1070] & -surfaces[1110]
lowerlevelwallcells[-1].name = "Auxiliary Wall high"

lowerlevelwallcells.append(openmc.Cell())
lowerlevelwallcells[-1].region = +surfaces[2010] & -surfaces[2020] & +surfaces[180] & -surfaces[200] & +surfaces[1050] & -surfaces[1060]
lowerlevelwallcells[-1].name = "Wing II Wall center"

lowerlevelwallcells.append(openmc.Cell())
lowerlevelwallcells[-1].region = +surfaces[2010] & -surfaces[2020] & +surfaces[190] & -surfaces[200] & +surfaces[1080] & -surfaces[1090]
lowerlevelwallcells[-1].name = "Wing I Wall auxiliary side"

lowerlevelwallcells.append(openmc.Cell())
lowerlevelwallcells[-1].region = +surfaces[2010] & -surfaces[2020] & +surfaces[190] & -surfaces[200] & +surfaces[1020] & -surfaces[1030]
lowerlevelwallcells[-1].name = "Wing I Wall weapon side"

lowerlevelwallcells.append(openmc.Cell())
lowerlevelwallcells[-1].region = +surfaces[2010] & -surfaces[2020] & +surfaces[200] & -surfaces[210] & +surfaces[1020] & -surfaces[1090]
lowerlevelwallcells[-1].name = "Wing I Wall high"

region_columnauxiliarylowouter = +surfaces[2020] & -surfaces[2060] & +surfaces[5] & -surfaces[10] & +surfaces[1070] & -surfaces[1080]
region_columnauxiliarylowinner = +surfaces[2020] & -surfaces[2060] & +surfaces[20] & -surfaces[30] & +surfaces[1070] & -surfaces[1080]
region_columnauxiliaryhighinner = +surfaces[2020] & -surfaces[2060] & +surfaces[180] & -surfaces[190] & +surfaces[1070] & -surfaces[1080]
region_columnauxiliaryhighouter = +surfaces[2020] & -surfaces[2060] & +surfaces[200] & -surfaces[210] & +surfaces[1070] & -surfaces[1080]

region_columnweaponlowouter = +surfaces[2020] & -surfaces[2060] & +surfaces[5] & -surfaces[10] & +surfaces[1030] & -surfaces[1040]
region_columnweaponlowinner = +surfaces[2020] & -surfaces[2060] & +surfaces[20] & -surfaces[30] & +surfaces[1030] & -surfaces[1040]
region_columnweaponhighinner = +surfaces[2020] & -surfaces[2060] & +surfaces[180] & -surfaces[190] & +surfaces[1030] & -surfaces[1040]
region_columnweaponhighouter = +surfaces[2020] & -surfaces[2060] & +surfaces[200] & -surfaces[210] & +surfaces[1030] & -surfaces[1040]

upperlevelwallcells = []

upperlevelwallcells.append(openmc.Cell())
upperlevelwallcells[-1].region = region_columnauxiliarylowouter
upperlevelwallcells[-1].name = "Column Auxiliary Low Outer"

upperlevelwallcells.append(openmc.Cell())
upperlevelwallcells[-1].region = region_columnauxiliarylowinner
upperlevelwallcells[-1].name = "Column Auxiliary Low Inner"

upperlevelwallcells.append(openmc.Cell())
upperlevelwallcells[-1].region = region_columnauxiliaryhighinner
upperlevelwallcells[-1].name = "Column Auxiliary High Inner"

upperlevelwallcells.append(openmc.Cell())
upperlevelwallcells[-1].region = region_columnauxiliaryhighouter
upperlevelwallcells[-1].name = "Column Auxiliary High Outer"

upperlevelwallcells.append(openmc.Cell())
upperlevelwallcells[-1].region = region_columnweaponlowouter
upperlevelwallcells[-1].name = "Column Weapon Low Outer"

upperlevelwallcells.append(openmc.Cell())
upperlevelwallcells[-1].region = region_columnweaponlowinner
upperlevelwallcells[-1].name = "Column Weapon Low Inner"

upperlevelwallcells.append(openmc.Cell())
upperlevelwallcells[-1].region = region_columnweaponhighinner
upperlevelwallcells[-1].name = "Column Weapon High Inner"

upperlevelwallcells.append(openmc.Cell())
upperlevelwallcells[-1].region = region_columnweaponhighouter
upperlevelwallcells[-1].name = "Column Weapon High Outer"

upperlevelwallcells.append(openmc.Cell())
upperlevelwallcells[-1].region = +surfaces[2020] & -surfaces[2060] & +surfaces[5] & -surfaces[210] & +surfaces[1020] & -surfaces[1030]
upperlevelwallcells[-1].name = "Wall Weapon Side"

upperlevelwallcells.append(openmc.Cell())
upperlevelwallcells[-1].region = +surfaces[2020] & -surfaces[2060] & +surfaces[5] & -surfaces[210] & +surfaces[1080] & -surfaces[1090]
upperlevelwallcells[-1].name = "Wall Auxiliary Side"

upperlevelroofcells = []
upperlevelroofcells.append(openmc.Cell())
upperlevelroofcells[-1].region = +surfaces[2020] & -surfaces[2040] & +surfaces[20] & -surfaces[190] & +surfaces[1000] & -surfaces[1020]
upperlevelroofcells[-1].name = "Weapon Roof"

upperlevelroofcells.append(openmc.Cell())
upperlevelroofcells[-1].region = +surfaces[2020] & -surfaces[2040] & +surfaces[20] & -surfaces[190] & +surfaces[1090] & -surfaces[1110]
upperlevelroofcells[-1].name = "Auxiliary Roof"

upperlevelroofcells.append(openmc.Cell())
upperlevelroofcells[-1].region = (+surfaces[2020] & -surfaces[2030] & +surfaces[5] & -surfaces[70] & +surfaces[1030] & -surfaces[1080]) & ~region_columnauxiliarylowouter & ~region_columnauxiliarylowinner & ~region_columnweaponlowouter & ~region_columnweaponlowinner
upperlevelroofcells[-1].name = "Floor low"

upperlevelroofcells.append(openmc.Cell())
upperlevelroofcells[-1].region = (+surfaces[2020] & -surfaces[2030] & +surfaces[140] & -surfaces[210] & +surfaces[1030] & -surfaces[1080]) & ~region_columnauxiliaryhighouter & ~region_columnauxiliaryhighinner & ~region_columnweaponhighouter & ~region_columnweaponhighinner
upperlevelroofcells[-1].name = "Floor high"

upperlevelroofcells.append(openmc.Cell())
upperlevelroofcells[-1].region = +surfaces[2060] & -surfaces[2070] & +surfaces[5] & -surfaces[210] & +surfaces[1020] & -surfaces[1090]
upperlevelroofcells[-1].name = "Upper Roof"

doors = []
doors.append(openmc.Cell())
doors[-1].region = +surfaces[2030] & -surfaces[2060] & +surfaces[5] & -surfaces[2] & +surfaces[1040] & -surfaces[1070]
doors[-1].name = "Main Door low"

doors.append(openmc.Cell())
doors[-1].region = +surfaces[2030] & -surfaces[2060] & +surfaces[28] & -surfaces[30] & +surfaces[1040] & -surfaces[1070]
doors[-1].name = "Inner Main Door low"

doors.append(openmc.Cell())
doors[-1].region = +surfaces[2030] & -surfaces[2060] & +surfaces[180] & -surfaces[182] & +surfaces[1040] & -surfaces[1070]
doors[-1].name = "Inner Main Door high"

doors.append(openmc.Cell())
doors[-1].region = +surfaces[2030] & -surfaces[2060] & +surfaces[208] & -surfaces[210] & +surfaces[1040] & -surfaces[1070]
doors[-1].name = "Main Door high"

doors.append(openmc.Cell())
doors[-1].region = +surfaces[2010] & -surfaces[2020] & +surfaces[37] &-surfaces[53] &  +surfaces[1030] & -surfaces[1032]
doors[-1].name = "Weapon Door low"

doors.append(openmc.Cell())
doors[-1].region = +surfaces[2010] & -surfaces[2020] & +surfaces[77] &-surfaces[93] &  +surfaces[1030] & -surfaces[1032]
doors[-1].name = "Weapon Door center low"

doors.append(openmc.Cell())
doors[-1].region = +surfaces[2010] & -surfaces[2020] & +surfaces[117] &-surfaces[133] & +surfaces[1030] & -surfaces[1032] 
doors[-1].name = "Weapon Door center high"

doors.append(openmc.Cell())
doors[-1].region = +surfaces[2010] & -surfaces[2020] & +surfaces[157] &-surfaces[173] & +surfaces[1030] & -surfaces[1032]
doors[-1].name = "Weapon Door high"

soilcells = []

soilcells.append(openmc.Cell())
soilcells[-1].region = +surfaces[2040] & +surfaces[2200] & +surfaces[20] & -surfaces[190] & +surfaces[1000] & -surfaces[1020]
soilcells[-1].name = "Weapon Area soil cover"

soilcells.append(openmc.Cell())
soilcells[-1].region = +surfaces[2040] & +surfaces[2210] & +surfaces[20] & -surfaces[190] & +surfaces[1090] & -surfaces[1110]
soilcells[-1].name = "Auxiliary Area soil cover"

soilcells.append(openmc.Cell())
soilcells[-1].region = +surfaces[2000] & +surfaces[2220] & +surfaces[5] & -surfaces[210] & -surfaces[1000]
soilcells[-1].name = "Weapon Area soil right dam"

soilcells.append(openmc.Cell())
soilcells[-1].region = +surfaces[2000] & -surfaces[2230] & +surfaces[5] & -surfaces[210] & +surfaces[1110]
soilcells[-1].name = "Auxiliary Area soil left dam"

soilcells.append(openmc.Cell())
soilcells[-1].region = +surfaces[2000] & +surfaces[2200] & +surfaces[5] & -surfaces[20] & +surfaces[1000] & -surfaces[1020]
soilcells[-1].name = "Weapon Area low dam"

soilcells.append(openmc.Cell())
soilcells[-1].region = +surfaces[2000] & +surfaces[2210] & +surfaces[5] & -surfaces[20] & +surfaces[1090] & -surfaces[1110]
soilcells[-1].name = "Auxiliary Area low dam"

soilcells.append(openmc.Cell())
soilcells[-1].region = +surfaces[2000] & +surfaces[2200] & +surfaces[190] & -surfaces[210] & +surfaces[1000] & -surfaces[1020]
soilcells[-1].name = "Weapon Area high dam"

soilcells.append(openmc.Cell())
soilcells[-1].region = +surfaces[2000] & +surfaces[2210] & +surfaces[190] & -surfaces[210] & +surfaces[1090] & -surfaces[1110]
soilcells[-1].name = "Auxiliary Area high dam"

soilcells.append(openmc.Cell())
soilcells[-1].region = +surfaces[2070] & -surfaces[2080] & +surfaces[5] & -surfaces[210] & +surfaces[1020] & -surfaces[1090]
soilcells[-1].name = "Center Area cover"


airCellBox = (+surfaces[2000] & -surfaces[2090] & +surfaces[1] & -surfaces[220] & +surfaces[990] & -surfaces[1120])
limitedAirCellBox = airCellBox
for cell in baseconcretecells + lowerlevelwallcells + upperlevelroofcells + upperlevelwallcells + doors + soilcells:
    limitedAirCellBox &= ~cell.region
    
airCell = openmc.Cell()
airCell.region = limitedAirCellBox
airCell.fill = airMat

for cell in baseconcretecells:
    cell.fill = wallMat

for cell in lowerlevelwallcells:
    cell.fill = wallMat

for cell in upperlevelwallcells:
    cell.fill = wallMat

for cell in upperlevelroofcells:
    cell.fill = wallMat

for cell in doors:
    cell.fill = steelMat

for cell in soilcells:
    cell.fill = soilMat

#*******************************************************************************
# fetter model

# Pit

cenSurface = openmc.Sphere(r = pitIR)
cenCell = openmc.Cell()
cenCell.region = -cenSurface
cenCell.name = "Center of weapon"

pitSurface = openmc.Sphere(r = pitOR)
pitCell = openmc.Cell()
pitCell.region = +cenSurface & -pitSurface
pitCell.name = "Pit"

# Reflector
refThickness = 2
refOR = pitOR + refThickness

refSurface = openmc.Sphere(r = refOR)
refCell = openmc.Cell()
refCell.region = +pitSurface & -refSurface
refCell.name = "Reflector"

# Tamper
tamThickness = 3
tamOR = refOR + tamThickness

tamSurface = openmc.Sphere(r = tamOR)
tamCell = openmc.Cell()
tamCell.region = +refSurface & -tamSurface
tamCell.name = "Tamper"

# Explosive
expThickness = 10
expOR = tamOR + expThickness

expSurface = openmc.Sphere(r = expOR)
expCell = openmc.Cell()
expCell.region = +tamSurface & -expSurface
expCell.name = "Conventional Explosive"

# Aluminum Case
aluThickness = 1
aluOR = expOR + aluThickness

aluSurface = openmc.Sphere(r = aluOR)
aluCell = openmc.Cell()
aluCell.region = +expSurface & -aluSurface
aluCell.name = "Aluminum Case"

weaponairCell = openmc.Cell()
weaponairCell.region = +aluSurface
weaponairCell.name = "Air around weapon"

pitCell.fill = wpuMat[weaponage] 
refCell.fill = berMat
tamCell.fill = tunMat
expCell.fill = hmxMat
aluCell.fill = aluMat
weaponairCell.fill = airMat


#*******************************************************************************
# cells & universes + source location

if fettersource:
    sourcex = surfaces[30].x0 + 50
    sourcey = surfaces[1010].y0 + 50
    sourcez = surfaces[2010].z0 + 50

    weaponsurface = openmc.Sphere(r = aluOR + 1, x0 = sourcex, y0 = sourcey, z0 = sourcez) # slightly bigger to avoid particles getting 'trapped' between universes, reaching MAX_EVENTS
    airCell.region &= +weaponsurface
else:
    sourcex = surfaces[30].x0 + 50
    sourcey = surfaces[1010].y0 + 50
    sourcez = surfaces[2010].z0 + 50

cells = baseconcretecells \
    + lowerlevelwallcells \
    + upperlevelwallcells \
    + upperlevelroofcells \
    + doors \
    + soilcells \
    + [airCell]
root = openmc.Universe(cells = cells)

if fettersource:
    fettercells = [cenCell, pitCell, refCell, tamCell, expCell, aluCell, weaponairCell]
    fetteruniverse = openmc.Universe(cells = fettercells)

    weaponcell = openmc.Cell(fill = fetteruniverse, region = -weaponsurface)
    weaponcell.translation = (sourcex, sourcey, sourcez)
    root.add_cell(weaponcell)

geometry = openmc.Geometry(root)
geometry.export_to_xml(os.path.join(basepath, "geometry.xml"))

# Plot
if plot:
    print("Plotting geometry - can take a few seconds")

    plt.figure(figsize=(8, 8))
    root.plot(origin = (sourcex, sourcey, sourcez),
              basis=('xz'),
              width=(300, 300),
              pixels=(600, 600),
              seed = 1)
    plt.title("Source View (xz)")
    plt.savefig(os.path.join(basepath, "source-view-xz.png"))
    plt.close()

    plt.figure(figsize=(8, 8))
    root.plot(origin = (sourcex, sourcey, sourcez),
              basis=('xy'),
              width=(300, 300),
              pixels=(600, 600),
              seed = 1)
    plt.title("Source View (xy)")
    plt.savefig(os.path.join(basepath, "source-view-xy.png"))
    plt.close()

    xfactor = 1.2
    yfactor = 1.2
    zfactor = 1.2

    plt.figure(figsize=(8, 8))
    root.plot(origin = (0 + xwidth / 2, 0 + ywidth / 2, thickA + 1),
              width=(xwidth * xfactor, ywidth * yfactor),
              seed = 1)
    plt.title("Lower Level floor plan")
    plt.savefig(os.path.join(basepath, "floor-plan-lower-level.png"))
    plt.close()

    plt.figure(figsize=(8, 8))
    root.plot(origin = (0 + xwidth / 2, 0 + ywidth / 2, thickA + distf + 1),
              width=(xwidth * xfactor, ywidth * yfactor),
              seed = 1)
    plt.title("Upper Level floor plan")
    plt.savefig(os.path.join(basepath, "floor-plan-upper-level.png"))
    plt.close()

    plt.figure(figsize=(8, 8))
    root.plot(origin = (0 + xwidth / 2, 0 + ywidth / 2, thickA + distf + thickC + 1),
              width=(xwidth * xfactor, ywidth * yfactor),
              seed = 1)
    plt.title("Upper Level doors floor plan")
    plt.savefig(os.path.join(basepath, "floor-plan-doors-upper-level.png"))
    plt.close()

    plt.figure(figsize=(8, 8))
    root.plot(origin = (0 + xwidth / 2, 0 + ywidth / 2, thickA + 1),
              width=(xwidth * xfactor, ywidth * yfactor),
              color_by='material')
    plt.title("Lower Level floor plan")
    plt.savefig(os.path.join(basepath, "floor-plan-lower-level-material.png"))
    plt.close()

    plt.figure(figsize=(8, 8))
    root.plot(origin = (0 + xwidth / 2, 0 + ywidth / 2, thickA + distf + 1),
              width=(xwidth * xfactor, ywidth * yfactor),
              color_by='material',
              seed = 1)
    plt.title("Upper Level floor plan")
    plt.savefig(os.path.join(basepath, "floor-plan-upper-level-material.png"))
    plt.close()

    plt.figure(figsize=(8, 8))
    root.plot(origin = (0 + xwidth / 2, 0 + ywidth / 2, thickA + distf + thickC + 1),
              width=(xwidth * xfactor, ywidth * yfactor),
              color_by='material',
              seed = 1)
    plt.title("Upper Level doors floor plan")
    plt.savefig(os.path.join(basepath, "floor-plan-doors-upper-level-material.png"))
    plt.close()

    yfactor = 1.5
    plt.figure(figsize=(8, 8))
    root.plot(origin = ((surfaces[20].x0 + surfaces[30].x0) / 2, (surfaces[1090].y0 + surfaces[1020].y0)/ 2, surfaces[2040].z0),
              basis = ('yz'),
              width=(ywidth * yfactor, zwidth * zfactor),
              color_by='cell',
              seed = 1)
    plt.title("Cut through soil")
    plt.savefig(os.path.join(basepath, "front-with-top-soil.png"))
    plt.close()

    yfactor = 1.5
    plt.figure(figsize=(8, 8))
    root.plot(origin = (surfaces[15].x0, (surfaces[1090].y0 + surfaces[1020].y0)/ 2, surfaces[2040].z0),
              basis = ('yz'),
              width=(ywidth * yfactor, zwidth * zfactor),
              color_by='cell',
              seed = 1)
    plt.title("Cut through soil")
    plt.savefig(os.path.join(basepath, "front-with-top-soil-low.png"))
    plt.close()

    yfactor = 1.5
    plt.figure(figsize=(8, 8))
    root.plot(origin = (surfaces[195].x0, (surfaces[1090].y0 + surfaces[1020].y0)/ 2, surfaces[2040].z0),
              basis = ('yz'),
              width=(ywidth * yfactor, zwidth * zfactor),
              color_by='cell',
              seed = 1)
    plt.title("Cut through soil")
    plt.savefig(os.path.join(basepath, "front-with-top-soil-high.png"))
    plt.close()
    
###############################################################################
# Source & Settings
###############################################################################

if fettersource:
    bounds = [sourcex - pitOR, sourcey - pitOR, sourcez - pitOR,
              sourcex + pitOR, sourcey + pitOR, sourcez + pitOR]
    uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable=True)
    subset = puagedf.loc[(slice(None), weaponage), :]
    normalizationconstant = sum(subset['wo'] * subset['sfneutrons'])
    source = []
    for iso in pudf.index:
        tmpsrc = openmc.Source(space = uniform_dist, particle = 'neutron')
        tmpsrc.energy = openmc.stats.Watt(a=1/pudf.loc[iso, 'watt-a'], b=pudf.loc[iso, 'watt-b'])
        tmpsrc.strength = puagedf.loc[(iso, age), 'sfneutrons'] * puagedf.loc[(iso, age), 'wo'] / normalizationconstant
        source.append(tmpsrc)
else:
    source = openmc.Source(space = openmc.stats.Point((sourcex, sourcey, sourcez)),
                           energy = openmc.stats.Discrete([2e6], [1]),
                           particle = 'neutron')

settings = openmc.Settings()
settings.run_mode = 'fixed source'
settings.inactive = 0
settings.batches = 10
settings.particles = particles
if survival:
    settings.survival_biasing = True
settings.source = source

settings.export_to_xml(os.path.join(basepath, "settings.xml"))

###############################################################################
# Tallies
###############################################################################

tallies = openmc.Tallies()

xmin = 0
xmax = xwidth
ymin = 0
ymax = ywidth
zmin = 0
zmax = zwidth

xstep = 100
tallyxmin = math.floor(xmin / xstep) * xstep
tallyxmax = math.ceil(xmax / xstep) * xstep 
ystep = 100
tallyymin = math.floor(ymin / ystep) * ystep
tallyymax = math.ceil(ymax / ystep) * ystep 
zstep = 100
tallyzmin = math.floor(zmin / zstep) * zstep
tallyzmax = math.ceil(zmax / zstep) * zstep 
xdim = int(abs(tallyxmax - tallyxmin) / xstep)
ydim = int(abs(tallyymax - tallyymin) / ystep)
zdim = int(abs(tallyzmax - tallyzmin) / zstep)

regularmesh = openmc.RegularMesh(name = "Regular Mesh")
regularmesh.dimension = (xdim, ydim, zdim)
regularmesh.lower_left = (tallyxmin, tallyymin, tallyzmin)
#regularmesh.upper_right = (tallyxmax, tallyymax, tallyzmax)
regularmesh.width = (xstep, ystep, zstep)

meshfilter = openmc.MeshFilter(regularmesh)
tally = openmc.Tally(name = "flux, regular mesh")
tally.filters = [meshfilter]
tally.scores = ["flux"]
tallies.append(tally)

meshsurfacefilter = openmc.MeshSurfaceFilter(regularmesh)
tally = openmc.Tally(name = "current, regular mesh")
tally.filters = [meshsurfacefilter]
tally.scores = ["current"]
tallies.append(tally)

tallies.export_to_xml(os.path.join(basepath, "tallies.xml"))
