import numpy as np
import openmc
import pandas as pd
import copy

from uncertainties import ufloat

def planeparameter3points(p1, p2, p3):
    if len(p1) != 3 and len(p2) != 3 and len(p3) != 3:
        print("All points need 3 dimensions")
        raise RuntimeError()
    p1 = np.array(p1)
    p2 = np.array(p2)
    p3 = np.array(p3)

    v1 = p3 - p1
    v2 = p2 - p1

    cp = np.cross(v1, v2)
    a, b, c = cp

    d = np.dot(cp, p3)

    return (a, b, c, d)

def checkcrosssections(cspath, materials):
    lib = openmc.data.DataLibrary.from_xml(cspath)
    for mat in materials:
        for nuclide in mat.nuclides:
            if lib.get_by_material(nuclide.name) is None:
                print("WARNING: Could not find {:s} in cross section library".format(nuclide.name))
                print("  The material contains {:e} {:s} of that nuclide.".format(nuclide.percent, nuclide.percent_type))
                print("  Will remove the nuclide from the material")
                mat.remove_nuclide(nuclide.name)

def createpuvecfrommat(material):
    puisotopes = ['Pu238', 'Pu239', 'Pu240', 'Pu241', 'Pu242', 'Am241']
    # totpumass = sum([material.get_mass_density(i) for i in puisotopes])
    totpumass = material.get_mass_density()
    puvec = puvector()
    for i in puisotopes:
        puvec.setwo(i, material.get_mass_density(i) / totpumass)
        print(i, material.get_mass_density(i) / totpumass)
    return puvec

def createwgpuvec():
    puvec = puvector()
    puvec.setwo('Pu238', 0.00005)
    puvec.setwo('Pu239', 0.993)
    puvec.setwo('Pu240', 0.060)
    puvec.setwo('Pu241', 0.0044)
    puvec.setwo('Pu242', 0.00015)
    puvec.setwo('Am241', 0)
    return puvec

class puvector:
    def __init__(self):
        pudict = {'name': ['Pu238', 'Pu239', 'Pu240', 'Pu241', 'Pu242', 'Am241'],
                  'sfneutrons-ufloat': [ufloat(2630, 140),
                                        ufloat(0.0152, 0.0029),
                                        ufloat(1031, 36),
                                        ufloat(0.0017233, 0.000035),
                                        ufloat(1722, 21),
                                        ufloat(1.47, 0.37)], # PhD Thesis Moritz
                  'watt-a': [1.17948e-6, 1.12963e-6, 1.25797e-6, 1.18698e-6, 1.22078e-6, 1.07179e-6], # Verbeke et al. 2010
                  'watt-b': [4.16933e-6, 3.80269e-6, 4.68927e-6, 4.15150e-6, 4.36668e-6, 3.46195e-6]} # Verbeke et al. 2010
        pudict['sfneutrons'] = [x.n for x in pudict['sfneutrons-ufloat']]
        pudict['mass'] = [openmc.data.atomic_mass(iso) for iso in pudict['name']]
        self.pudf = pd.DataFrame(pudict)
        self.pudf['wo'] = 0
        self.pudf.set_index('name', inplace = True)

    def setwo(self, iso, fraction):
        self.pudf.loc[iso, 'wo'] = fraction

    def createagedvector(self, ages = [0, 5, 10, 15, 20, 25, 30, 35, 40]):
        self.puagedf = pd.DataFrame(columns = ['name', 'age', 'wo', 'sfneutrons', 'watt-a', 'watt-b'])
        self.puagedf.set_index(['name', 'age'], inplace = True)

        # all times in years
        lbda = np.log(2)/14.329 # T_1/2 from PhD Thesis Moritz
        totmol = sum(self.pudf['wo'] / self.pudf['mass']) + 0.002 / openmc.data.atomic_weight('O')

        Npu0 = self.pudf.loc['Pu241', 'wo'] / openmc.data.atomic_mass('Pu241') / totmol 
        for age in ages:
            agelist = [age] * 6
            Npu = Npu0 * np.exp(-lbda * age)
            Nam = Npu0 * (1 - np.exp(-lbda * age))
            puadict = self.pudf.to_dict()
            puadict['wo']['Pu241'] = Npu * totmol * openmc.data.atomic_mass('Pu241')
            puadict['wo']['Am241'] = Nam * totmol * openmc.data.atomic_mass('Am241')
            puadict['age'] = {iso: age for iso in puadict['wo']}

            tempdf = pd.DataFrame(puadict)
            tempdf.index = tempdf.index.set_names(['name'])    
            tempdf.reset_index(inplace=True)
            tempdf.set_index(['name', 'age'], inplace=True)
            self.puagedf = pd.concat((self.puagedf, tempdf))

        # self.puagedf.set_index(['name', 'age'], inplace=True)
        
