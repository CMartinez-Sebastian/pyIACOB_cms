#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: C. Martinez-Sebastian
"""
import pandas as pd
import os
import pickle
import tarfile
import numpy as np
import time
from scipy.interpolate import interp1d


class Microturbulence:
    
    def __init__(self, ew, wvl, flx):
        '''
        Parameters
        ----------
        ew : float,
            Equivalent Width.
        wvl : array,
            Wavelenth of a given profile line
        flx : array,
            Flux of a given profile line.

        Returns
        -------
        None.
        '''
        self.EW = ew
        self.wvl = wvl
        self.flx = flx

class Abundance:
    def __init__(self):
        self.microturbulences = {}
    
    def add_microturbulence(self, microturbulence, name):
        '''
        Parameters
        ----------
        microturbulence : class,
            Class Microturbulence
        name : str,
            Name of the corresponding microturbulence (normally in format VT000).

        Returns
        -------
        None.
        '''
        self.microturbulences[name] = microturbulence
        
    def get_microturbulence(self, name):
        '''
        Parameters
        ----------
        name : str,
            Name of the desired microturbulence (normally in format VT000)

        Returns
        -------
        Microturbulence : class
        '''
        return self.microturbulences[name]

class Line:
    def __init__(self):
        self.abundances = {}

    def add_abundance(self, abundance, name):
        '''
        Parameters
        ----------
        abundance : Class,
            Class Abundance.
        name : str,
            Name of the corresponding microturbulence (normally in decimal-negative format)

        Returns
        -------
        None
        '''
        self.abundances[name] = abundance
    def get_abundance(self, name):
        '''
        Parameters
        ----------
        name : 
            Name of the desired abundance.

        Returns 
        -------
        Abundance : Class,
        '''
        return self.abundances[name]
    
    def get_microturbulence(self, name_abundance, name_microturbulence):
        '''
        Parameters
        ----------
        name_abundance : str,
            Name of the desired abundance
        name_microturbulence : str,
            Name of the desired microturbulence for the given abundance

        Returns
        -------
        Microturbulence : class, 
        '''
        abundance = self.get_abundance(name_abundance)
        return abundance.microturbulences(name_microturbulence)
    
class Header:
    def __init__(self, hd):
        '''
        Parameters
        ----------
        hd : list,
            List with the values from INDAT file (except abundance)

        Returns
        -------
        None.

        '''
        self.name = hd[0]
        self.T_eff, self.log_g, self.R_Rsun = hd[1],hd[2],hd[3]
        self.R_max, self.T_min = hd[4],hd[5]
        self.M_dot, self.v_inf, self.beta = hd[6],hd[7],hd[8]
        self.Y_He = hd[9]
        self.v_turb, self.Z_t = hd[10],hd[11]
        self.cl1,self.cl2,self.cl3 = hd[12],hd[13],hd[14]
    def add_abundance(self,abundance):
        '''
        Possibility of adding the abundance
        
        Parameters
        ----------
        abundance : float,
            Abundance value

        Returns
        -------
        None.

        '''
        self.abundance = abundance

class Model:
    def __init__(self):
        self.lines = {}
        self.header = None
        
    def add_line(self, line, name):
        '''
        Adds a given line to the model 
        
        Parameters
        ----------
        line : class,
            Class Line
        name : str,
            Line name.

        Returns
        -------
        None.
        '''
        self.lines[name] = line
    
    def add_header(self,header):
        '''
        Adds a header to the model

        Parameters
        ----------
        header : list,
            Parameter of Header class.

        Returns
        -------
        None.
        '''
        self.header = header
        
    def get_line(self, name):
        return self.lines[name]

    def get_abundance(self, name_line, name_abundance):
        abundance = self.get_line(name_line)
        return abundance.get_abundance(name_abundance)
    
    def get_profile(self, name_line, name_abundance, name_micr):
        '''
        Parameters
        ----------
        name_line : str,
            Name of a given line.
        name_abundance : str,
            Name of the desired abundance (normally in negative-decimal format)
        name_micr : str,
            Name of the desired microturbulence (normally in VT000 format).

        Returns
        -------
        Microturbulence : class,
           Corresponding profile for the given line, with the given abundance, and 
           microturbulence
        '''
        abudance = self.get_abundance(name_line, name_abundance)
        return abudance.get_microturbulence(name_micr)

def OUT_to_pkl_lines(model, elem): # To save the original lines in a different format.
    '''
    Function to convert from the original OUT format from FASTWIND (ascii) to 
    pickle python format (binary).

    It saves a pkl file in "original/'model+elem'_lines.pkl", with the structure
    of class Model.
    
    Parameters
    ----------
    model: str,
        String with the name of a given model. It assumes a name of the form
        T400g390He10Q130b15CNO (INCLUDING THE FINAL CNO, WHICH IS REMOVED IN 
        THE PROGRAM)
        
    elem: str,
        Name of the element to convert. It is used to give the saving name and
        to look for the abundance and line files: abundance_'elem'.dat and
        lines_'elem'.dat
        
    Warnings
    --------
        You have to make attention to the models name and the folders in which 
        the different files are looked for: lines, abundances and microturbulences
        files are looked for in the same folder in which we are working (supposed
        to be in the original grid folder) or in 'grid_components' folder (previous
        level).
    
    Returns
    -------
    None
    '''
    
    mod = Model()
    ########################## OJO! HABRA QUE LEERLOS DE UN DAT ##########################
    try: 
        metallicities = pd.read_csv('abundances_'+elem+'.dat', header=None,delim_whitespace=True)
        microturbulences = pd.read_csv('microturbulences.dat', header=None,delim_whitespace=True)
        LINES = pd.read_csv('lines_'+elem+'.dat', header=None,delim_whitespace=True)
    except: 
        metallicities = pd.read_csv('../grid_components/abundances_'+elem+'.dat', header=None,delim_whitespace=True)
        microturbulences = pd.read_csv('../grid_components/microturbulences.dat', header=None,delim_whitespace=True)
        LINES = pd.read_csv('../grid_components/lines_'+elem+'.dat', header=None,delim_whitespace=True)
    ########################## OJO! HABRA QUE LEERLOS DE UN DAT ##########################
    possible_models = os.listdir('.')
    start = time.time()
    for line in LINES[0]:
        lin = Line()
        for ab in metallicities[1]:  # Reading all OUTs with the different given metallicities/abundancies. 
            abundance = Abundance()
            lin.add_abundance(abundance,str(ab))
            del abundance
        mod.add_line(lin, line)
    for idx in metallicities.index:
        Z, ab = metallicities[0][idx], metallicities[1][idx]
        if model+Z in possible_models:
            for line in LINES[0]:  # Reading all OUTs with the different given metallicities/abundancies. 
               for VT in microturbulences[0]:
                   try:                        # Working with compressed files
                       A = model+Z+'/OUT.'+line+'_'+VT        
                       DATA = pd.read_csv(A, header=None,delim_whitespace=True,usecols=[0,2,4],verbose=False,names=['EW','wvl','flx'])
                       wvl,flx = pd.to_numeric(DATA.wvl,errors='coerce'),pd.to_numeric(DATA.flx,errors='coerce')
                       wvl,flx = wvl[~np.isnan(wvl)],flx[~np.isnan(flx)] 
                       EW = float(DATA.EW[len(flx)])
                       del DATA
                       micro = Microturbulence(EW*-1e-3,wvl,flx) # Saving the different microturbulences
                       mod.lines[line].abundances[str(ab)].add_microturbulence(micro, VT)
                       del micro,EW,wvl,flx
                   except: 
                       micro = Microturbulence(None,None,None) # Saving the different microturbulences
                       mod.lines[line].abundances[str(ab)].add_microturbulence(micro, VT)
                       del micro
               if not mod.header:
                   A = open(model+Z+'/INDAT.DAT','r')
                   lines = A.readlines()
                   A.close()
                   model_name = lines[0].split()[0][1:]
                   T_eff, log_g, R_Rsun = float(lines[3].split()[0]),float(lines[5].split()[1]),float(lines[5].split()[2])
                   R_max, T_min = float(lines[4].split()[0]),float(lines[5].split()[1])
                   M_dot, v_inf, beta = float(lines[5].split()[0]),float(lines[5].split()[2]),float(lines[5].split()[3])
                   Y_He = float(lines[6].split()[0])
                   v_turb, Z_t = float(lines[8].split()[0]),float(lines[8].split()[1])
                   cl1,cl2,cl3 = float(lines[10].split()[0]),float(lines[10].split()[1]),float(lines[10].split()[2])
                   info = [model_name,T_eff,log_g,R_Rsun,R_max,T_min,M_dot,v_inf,beta,Y_He,v_turb,Z_t,cl1,cl2,cl3]
                   header = Header(info)
                   mod.add_header(header)
        else:
            try:
                tar = tarfile.open(model+Z+'.tgz','r')
                for line in LINES[0]:  # Reading all OUTs with the different given metallicities/abundancies. 
                   for VT in microturbulences[0]:
                       try:                        # Working with compressed files
                           file = model+Z+'/OUT.'+line+'_'+VT        
                           A =  tar.extractfile(file)
                           DATA = pd.read_csv(A, header=None,delim_whitespace=True,usecols=[0,2,4],verbose=False,names=['EW','wvl','flx'])
                           wvl,flx = pd.to_numeric(DATA.wvl,errors='coerce'),pd.to_numeric(DATA.flx,errors='coerce')
                           wvl,flx = wvl[~np.isnan(wvl)],flx[~np.isnan(flx)] 
                           EW = float(DATA.EW[len(flx)])
                           del DATA
                           micro = Microturbulence(EW*-1e-3,wvl,flx) # Saving the different microturbulences
                           mod.lines[line].abundances[str(ab)].add_microturbulence(micro, VT)
                           del micro,EW,wvl,flx,file
                       except: 
                           micro = Microturbulence(None,None,None) # Saving the different microturbulences
                           mod.lines[line].abundances[str(ab)].add_microturbulence(micro, VT)
                           del micro
                   if not mod.header:
                       A = tar.extractfile(model+Z+'/INDAT.DAT')
                       lines = A.readlines()
                       A.close()
                       model_name = lines[0].split()[0][1:]
                       T_eff, log_g, R_Rsun = float(lines[3].split()[0]),float(lines[5].split()[1]),float(lines[5].split()[2])
                       R_max, T_min = float(lines[4].split()[0]),float(lines[5].split()[1])
                       M_dot, v_inf, beta = float(lines[5].split()[0]),float(lines[5].split()[2]),float(lines[5].split()[3])
                       Y_He = float(lines[6].split()[0])
                       v_turb, Z_t = float(lines[8].split()[0]),float(lines[8].split()[1])
                       cl1,cl2,cl3 = float(lines[10].split()[0]),float(lines[10].split()[1]),float(lines[10].split()[2])
                       info = [model_name,T_eff,log_g,R_Rsun,R_max,T_min,M_dot,v_inf,beta,Y_He,v_turb,Z_t,cl1,cl2,cl3]
                       header = Header(info)
                       mod.add_header(header)
                       del header
                tar.close()  
            except: print('Not '+model+Z)
    end = time.time()
    print('Filling class time: ',end-start)
    #with open('../iac_iacobgrid_HHeCNO_Z100_He10_v106_pkl/original/'+model[:-3]+elem\
    #          +'_lines.pkl', 'wb') as outp:
    #    pickle.dump(mod, outp, pickle.HIGHEST_PROTOCOL)
    with open('../iac_iacobgrid_HHeCNO_Z100_He10_v106_pkl/original/'+model[:-3]+elem\
              +'_lines.pkl', 'wb') as outp:
        pickle.dump(mod, outp, pickle.HIGHEST_PROTOCOL)
        
    outp.close()   
    del mod
    return None#mod
print('class_FASTWIND read')

'''

with open('/media/cms/Elements/PhD/DATA/FASTWIND/iac_gbatgrid_HHeCNO_Z100/iac_gbatgrid_HHeCNO_Z100_He10_pkl/original/T550g400He10Q1230b15C_lines.pkl','rb') as inp:
    C = pickle.load(inp, encoding='latin1')


import matplotlib.pyplot as plt
plt.plot(C.lines['CIV580111'].abundances['-3.55'].microturbulences['VT009'].wvl,\
         C.lines['CIV580111'].abundances['-3.55'].microturbulences['VT009'].flx)

    
    
#%%


#with open(folder_original_pkl+'T580g410He10Q121b15C_lines.pkl', 'rb') as inp:
#    ex = pickle.load(inp, encoding='latin1')

M = 'T600g450He10Q140b15CNO' 
try: os.chdir("/net/nas/proyectos/hots/AAA_WORKING_2030/grids/fastwind/FW10.6/iac_iacobgrid_HHeCNO_Z100/iac_iacobgrid_HHeCNO_Z100_He10_v106")
except: os.chdir("/media/cms/Elements/PhD/DATA/FASTWIND/iac_iacobgrid_HHeCNO_Z100_He10_v106")
elems = ['C','N','O','H_He']
for el in elems:
    start_fin = time.time()

    globals()['EXAMPLE_'+el] = OUT_to_pkl_lines(M, el)
    end_fin = time.time()
    print(el,end_fin - start_fin)

# In[proof]

try: os.chdir("/net/nas/proyectos/hots/AAA_WORKING_2030/grids/fastwind/FW10.6/iac_iacobgrid_HHeCNO_Z100/iac_iacobgrid_HHeCNO_Z100_He10_v106")
except: os.chdir("/media/cms/Elements/PhD/DATA/FASTWIND/iac_iacobgrid_HHeCNO_Z100_He10_v106")
elems = ['C','N','O','H_He']
for el in elems:
    start_fin = time.time()

    globals()['EXAMPLE_'+el] = OUT_to_pkl_lines('T600g450He10Q140b15CNO', el)
    end_fin = time.time()
    print(el,end_fin - start_fin)

#%%

wvl1 = EXAMPLE_C.lines['CIV580111'].abundances['-2.2'].microturbulences['VT007'].wvl
wvl2 = EXAMPLE_C.lines['CIV580111'].abundances['-3.7'].microturbulences['VT007'].wvl
wvl3 = EXAMPLE_C.lines['CIV580111'].abundances['-4.3'].microturbulences['VT007'].wvl

flx1 = EXAMPLE_C.lines['CIV580111'].abundances['-2.2'].microturbulences['VT007'].flx
flx2 = EXAMPLE_C.lines['CIV580111'].abundances['-3.7'].microturbulences['VT007'].flx
flx3 = EXAMPLE_C.lines['CIV580111'].abundances['-4.3'].microturbulences['VT007'].flx

import matplotlib.pyplot as plt
plt.plot(wvl1,flx1)
plt.plot(wvl2,flx2)
plt.plot(wvl3,flx3)


DATAp000 = pd.read_csv('T600g450He10Q140b15CNOp000/OUT.CIV580111_VT007',\
                       header=None,delim_whitespace=True,usecols=[2,4],verbose=False,names=['wvl','flx'])

DATAp150 = pd.read_csv('T600g450He10Q140b15CNOp150/OUT.CIV580111_VT007',\
                       header=None,delim_whitespace=True,usecols=[2,4],verbose=False,names=['wvl','flx'])

DATAm060 = pd.read_csv('T600g450He10Q140b15CNOm060/OUT.CIV580111_VT007',\
                       header=None,delim_whitespace=True,usecols=[2,4],verbose=False,names=['wvl','flx'])

plt.plot(pd.to_numeric(DATAp150.wvl,errors='coerce'),pd.to_numeric(DATAp150.flx,errors='coerce'))
plt.plot(pd.to_numeric(DATAp000.wvl,errors='coerce'),pd.to_numeric(DATAp000.flx,errors='coerce'))
plt.plot(pd.to_numeric(DATAm060.wvl,errors='coerce'),pd.to_numeric(DATAm060.flx,errors='coerce'))
'''