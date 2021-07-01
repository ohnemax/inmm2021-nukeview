import openmc
import json
import os
import matplotlib.pyplot as plt
import math

basepath = "sost-20210626-100k"
particles = 100000

cspath = "/openmc/openmc-data/v0.12/lib80x_hdf5/cross_sections.xml"

if not os.path.exists(basepath):
    os.mkdir(basepath)

# dump basic settings to json file
with open(os.path.join(basepath, "calculation.json"), 'w') as f:
    json.dump({'particles': particles,
               'basepath': basepath,
               'cspath': cspath},
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

materiallist = [steelMat, concreteRegularMat, wallMat, airMat]
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
                     0: distp,
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
surfaces[2] = openmc.XPlane(surfaces[0].x0 + thickD, surface_id = 2)
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

# x surfaces for lower level (auxiliary)

# y surfaces for lower level (main)
surfaceincrements = {990: 0,
                     1000: distp,
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


# cells: base concrete
baseconcretecells = []
baseconcretecells.append(openmc.Cell())
baseconcretecells[-1].region = +surfaces[2000] & -surfaces[2010] & +surfaces[0] & -surfaces[20] & +surfaces[1020] & -surfaces[1090]
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
lowerlevelwallcells[-1].region = +surfaces[2010] & -surfaces[2020] & +surfaces[0] & -surfaces[10] & +surfaces[1020] & -surfaces[1090]
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

region_columnauxiliarylowouter = +surfaces[2020] & -surfaces[2060] & +surfaces[0] & -surfaces[10] & +surfaces[1070] & -surfaces[1080]
region_columnauxiliarylowinner = +surfaces[2020] & -surfaces[2060] & +surfaces[20] & -surfaces[30] & +surfaces[1070] & -surfaces[1080]
region_columnauxiliaryhighinner = +surfaces[2020] & -surfaces[2060] & +surfaces[180] & -surfaces[190] & +surfaces[1070] & -surfaces[1080]
region_columnauxiliaryhighouter = +surfaces[2020] & -surfaces[2060] & +surfaces[200] & -surfaces[210] & +surfaces[1070] & -surfaces[1080]

region_columnweaponlowouter = +surfaces[2020] & -surfaces[2060] & +surfaces[0] & -surfaces[10] & +surfaces[1030] & -surfaces[1040]
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
upperlevelwallcells[-1].region = +surfaces[2020] & -surfaces[2060] & +surfaces[0] & -surfaces[210] & +surfaces[1020] & -surfaces[1030]
upperlevelwallcells[-1].name = "Wall Weapon Side"

upperlevelwallcells.append(openmc.Cell())
upperlevelwallcells[-1].region = +surfaces[2020] & -surfaces[2060] & +surfaces[0] & -surfaces[210] & +surfaces[1080] & -surfaces[1090]
upperlevelwallcells[-1].name = "Wall Auxiliary Side"

upperlevelroofcells = []
upperlevelroofcells.append(openmc.Cell())
upperlevelroofcells[-1].region = +surfaces[2020] & -surfaces[2040] & +surfaces[20] & -surfaces[190] & +surfaces[1000] & -surfaces[1020]
upperlevelroofcells[-1].name = "Weapon Roof"

upperlevelroofcells.append(openmc.Cell())
upperlevelroofcells[-1].region = +surfaces[2020] & -surfaces[2040] & +surfaces[20] & -surfaces[190] & +surfaces[1090] & -surfaces[1110]
upperlevelroofcells[-1].name = "Auxiliary Roof"

upperlevelroofcells.append(openmc.Cell())
upperlevelroofcells[-1].region = (+surfaces[2020] & -surfaces[2030] & +surfaces[0] & -surfaces[70] & +surfaces[1030] & -surfaces[1080]) & ~region_columnauxiliarylowouter & ~region_columnauxiliarylowinner & ~region_columnweaponlowouter & ~region_columnweaponlowinner
upperlevelroofcells[-1].name = "Floor low"

upperlevelroofcells.append(openmc.Cell())
upperlevelroofcells[-1].region = (+surfaces[2020] & -surfaces[2030] & +surfaces[140] & -surfaces[210] & +surfaces[1030] & -surfaces[1080]) & ~region_columnauxiliaryhighouter & ~region_columnauxiliaryhighinner & ~region_columnweaponhighouter & ~region_columnweaponhighinner
upperlevelroofcells[-1].name = "Floor high"

upperlevelroofcells.append(openmc.Cell())
upperlevelroofcells[-1].region = +surfaces[2060] & -surfaces[2070] & +surfaces[0] & -surfaces[210] & +surfaces[1020] & -surfaces[1090]
upperlevelroofcells[-1].name = "Upper Roof"

doors = []
doors.append(openmc.Cell())
doors[-1].region = +surfaces[2030] & -surfaces[2060] & +surfaces[0] & -surfaces[2] & +surfaces[1040] & -surfaces[1070]
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


airCellBox = (+surfaces[2000] & -surfaces[2090] & +surfaces[1] & -surfaces[220] & +surfaces[990] & -surfaces[1120])
limitedAirCellBox = airCellBox
for cell in baseconcretecells + lowerlevelwallcells + upperlevelroofcells + upperlevelwallcells + doors:
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
    
cells = baseconcretecells \
    + lowerlevelwallcells \
    + upperlevelwallcells \
    + upperlevelroofcells \
    + doors \
    + [airCell]
    
root = openmc.Universe(cells = cells)
geometry = openmc.Geometry(root)
geometry.export_to_xml(os.path.join(basepath, "geometry.xml"))

# Plot
# print("Plotting geometry - can take a few seconds")
# xfactor = 1.2
# yfactor = 1.2
# plt.figure(figsize=(8, 8))
# root.plot(origin = (0 + xwidth / 2, 0 + ywidth / 2, thickA + 1), width=(xwidth * xfactor, ywidth * yfactor))
# plt.title("Lower Level floor plan")
# plt.savefig(os.path.join(basepath, "floor-plan-lower-level.png"))

# plt.figure(figsize=(8, 8))
# root.plot(origin = (0 + xwidth / 2, 0 + ywidth / 2, thickA + distf + 1), width=(xwidth * xfactor, ywidth * yfactor))
# plt.title("Upper Level floor plan")
# plt.savefig(os.path.join(basepath, "floor-plan-upper-level.png"))

# plt.figure(figsize=(8, 8))
# root.plot(origin = (0 + xwidth / 2, 0 + ywidth / 2, thickA + distf + thickC + 1), width=(xwidth * xfactor, ywidth * yfactor))
# plt.title("Upper Level doors floor plan")
# plt.savefig(os.path.join(basepath, "floor-plan-doors-upper-level.png"))

# plt.figure(figsize=(8, 8))
# root.plot(origin = (0 + xwidth / 2, 0 + ywidth / 2, thickA + 1),
#           width=(xwidth * xfactor, ywidth * yfactor),
#           color_by='material')
# plt.title("Lower Level floor plan")
# plt.savefig(os.path.join(basepath, "floor-plan-lower-level-material.png"))

# plt.figure(figsize=(8, 8))
# root.plot(origin = (0 + xwidth / 2, 0 + ywidth / 2, thickA + distf + 1),
#           width=(xwidth * xfactor, ywidth * yfactor),
#           color_by='material')
# plt.title("Upper Level floor plan")
# plt.savefig(os.path.join(basepath, "floor-plan-upper-level-material.png"))

# plt.figure(figsize=(8, 8))
# root.plot(origin = (0 + xwidth / 2, 0 + ywidth / 2, thickA + distf + thickC + 1),
#           width=(xwidth * xfactor, ywidth * yfactor),
#           color_by='material')
# plt.title("Upper Level doors floor plan")
# plt.savefig(os.path.join(basepath, "floor-plan-doors-upper-level-material.png"))

###############################################################################
# Source & Settings
###############################################################################

sourcex = surfaces[30].x0 + 50
sourcey = surfaces[1010].y0 + 50
sourcez = surfaces[2010].z0 + 50

source = openmc.Source(space = openmc.stats.Point((sourcex, sourcey, sourcez)),
                       energy = openmc.stats.Discrete([2e6], [1]),
                       particle = 'neutron')


settings = openmc.Settings()
settings.run_mode = 'fixed source'
settings.inactive = 0
settings.batches = 100
settings.particles = particles
settings.source = source

settings.export_to_xml(os.path.join(basepath, "settings.xml"))

###############################################################################
# Tallies
###############################################################################

tallies = openmc.Tallies()

xstep = 100
xmax = math.ceil(xwidth / 100) * 100 + 1
ystep = 100
ymax = math.ceil(ywidth / 100) * 100 + 1
zstep = 100
zmax = math.ceil(zwidth / 100) * 100 + 1

mainmesh = openmc.RectilinearMesh(name = "Main Mesh")
mainmesh.x_grid = range(0, xmax, xstep)
mainmesh.y_grid = range(0, ymax, ystep)
mainmesh.z_grid = range(0, zmax, zstep)

meshfilter = openmc.MeshFilter(mainmesh)
tally = openmc.Tally(name = "flux")
tally.filters = [meshfilter]
tally.scores = ["flux"]
tallies.append(tally)

tallies.export_to_xml(os.path.join(basepath, "tallies.xml"))
