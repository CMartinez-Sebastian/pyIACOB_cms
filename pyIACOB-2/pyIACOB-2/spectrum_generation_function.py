#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: C. Martinez-Sebastian

Different useful functions to generate synthetic spectra
"""

import pandas as pd
import numpy as np

import pickle
from class_FASTWIND import *
from math import fsum
from scipy.interpolate import interp1d
import scipy.constants as cte
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator,LinearLocator
import os
from matplotlib.backends.backend_pdf import PdfPages

import sys
sys.path.append('/home/cms/Desktop/PhD/PROJECT/pyIACOB-2')
import line_broadening as lb
from db import dir

mod = 'new'
print(mod)


def find_nearest(array, value):
    '''
    Parameters
    ----------
    array : array-like,
        Array in which we want to find the nearest value to a given one.
    value : float,
        Value for which we want to find the nearest value.

    Returns
    -------
    Nearest value IN the array TO the value

    '''
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def spectrum_generation(model,elems,elems_abundance,VT,wv1=3800,wv2=7000,stp=0.01,info=False,\
                        folder_original_pkl=dir().fastwinddir+dir().fastwinddir.split('/')[-2]+'_pkl/original/'):
    '''
    Function to create the spectrum from a given element with the considered lines. 
    It is NOT broadened. 
    
    Parameters
    ----------
    model : str,
        Model name of the form 'T360g390He10Q130b10CNO' (INCLUDING the CNO specification, which is the only computed at this time)
        
    elems : str,
        Name of the elements to be added to the model name (after removing from it CNO). At the moment: C, N, O, H_He
        
    elems_abundance : dictionary. 
        Element abundance. It has to be in the pkl. 
        Accept values in str or float. 
        If values are in format log[X/H]+12, it substracts 12
        
    VT : str or float/int,
        Microturbulence to be considered. It has to be in the pkl, and usually has the form 'VT007'.
        If type is not str, it looks for the closest value and write the corresponding format.
        
    wv1 : float, optional
        First wavelength to be considered. The default is 3800.
        
    wv2 : float, optional
        Last wavelength to be considered. The default is 7000.
        
    stp : float, optional
        step size for the wavelength base. The default is 0.01.
        
    info : Bool, optional 
        If true, it prints the line, maximum line flux, and the wavelength corresponding
        to that. Default is False
        
    folder_original_pkl : str, optional
        Folder to find the pkl files from the different elements and models. 
        The default is fastwinddir+*fastwind_grid*_pkl/original/'.

    Returns
    -------
    wave_range : 
        Wavelength array 
        
    flux_range
        Flux array (normalized to 1)
    '''
    os.chdir(dir().fastwinddir)
    with open('non_repeated_lines_Cnew.txt') as f:
        lines = f.readlines()
    f.close()

    considered_lines = []
    for l in lines:
        try: 
            if int(l.split()[1].replace(",",""))==1:
                considered_lines.append(l.split()[0].replace(",",""))
        except: pass
    
    ########################## Interpolated wave base #########################
    wave_range = np.around(np.arange(wv1,wv2,stp),2)
    flux_range = np.zeros(len(wave_range))
    if type(VT) != str:
        print(VT)
        try:
            microturbulences = pd.read_csv(dir().modeldir+'grid_components/microturbulences.dat', header=None,delim_whitespace=True)      
            FW_mic = [float(M[-3:]) for M in microturbulences[0]]
            VT = find_nearest(FW_mic,VT)
            print('Microturbulence is',VT)
            VT = 'VT000'[:-len(str(int(np.rint(VT))))]+str(int(np.rint(VT)))
        except: 
            print('You should use a suitable microturbulence. Not found the correct VT')
            return None
    for elem in elems:
        ########################## Reading pkl files ##########################
        with open(folder_original_pkl + model[:-3] + elem + '_lines.pkl', 'rb') as inp:
            profiles = pickle.load(inp, encoding='latin1')
            inp.close()
        ############ Getting model lines parameters for the element ###########
        LINES = pd.read_csv(dir().modeldir+'grid_components/lines_'+elem+'.dat', header=None,delim_whitespace=True)
        metallicity = elems_abundance[elem]
        if float(metallicity)>0: 
            metallicity -= 12
            metallicity = round(metallicity,2)
            metallicities = pd.read_csv(dir().modeldir+'grid_components/abundances_'+elem+'.dat', header=None,delim_whitespace=True)[1] # New????
            metallicity = find_nearest(metallicities,metallicity)
        ##### Interpolation of the different lines to the same base & sum #####
        for line in LINES[0]:
            if line in considered_lines:
                #print(line,str(metallicity),VT)
                line_wavelength = profiles.get_profile(line,str(metallicity),VT).wvl 
                if any(line_wavelength>=wv1) and any(line_wavelength<=wv2):      
                    line_flux = profiles.get_profile(line,str(metallicity),VT).flx 
                    if info is True:
                        print(line,max(line_flux),line_wavelength[np.argmax(line_flux)])
                    f1 = interp1d(line_wavelength,line_flux-1, kind='linear')  
                    mn, Mx = np.nanmin(line_wavelength),np.nanmax(line_wavelength)
                    mask = (wave_range>=mn) & (wave_range<=Mx)
                    masked_wave = wave_range[mask]
                    flx = f1(masked_wave)
                    flux_range[mask] += flx
    return wave_range,flux_range+1

 
class synthetic_model():
    def __init__(self, E_C, E_N, E_O, vmic, model, R=85000, vsini=20, vmac=5, beta=1.5, **kwargs):
        '''
        # def synthetic_model(star, E_C, E_N, E_O, vmic, model, R = 85000, vsini = 20, vmac = 5, beta = 1.5, **kwargs):
        Function to generate the complete synthetic spectra (broadened, etc)
    
        Parameters
        ----------
        ### star : str            ### Removed from the new version
        ###    Name of the star
            
        E_C, E_N, E_O : float
            Abundances of the different elements. For He, it uses E_N. 
            If values are in format log[X/H]+12, it substracts 12
            It has the form log(X/H) + 12
            
        vmic : str or float/int,
            Microturbulence to be considered. It has to be in the pkl, and usually has the form 'VT007'.
            If type is not str, it looks for the closest value and write the corresponding format.     
            
        model : str,
            Name of the model to consider. Has to include (end by) 'CNO' extension   
        
        R : float or int, optional
            Resolution to be interpolated to. Default is 85000
            
        vsin : float, optional
            Projected rotational velocity in km/s. The default is None.
    
            
        vmac : float, optional    
            Macroturbulence in km/s. The default is None.
            
        beta : float, optional
            Darkening parameter in the rotational broadening. It has to be between 0 and 1.5. The default is 1.5.
    
        **Kwargs :
            Arguments for function "spectrum_generation"
    
        Returns
        -------
        Wave : array-like
            Wavelength of the generated spectra
            
        Flux : array-like
            Flux of the generated spectra
    
        '''
        E_H_He = E_N
        ##### Reading possible microturbulence and finding the closest to the given #####
        if type(vmic) != str:
            try:
                microturbulences = pd.read_csv(dir().modeldir+'grid_components/microturbulences.dat', header=None,delim_whitespace=True)      
                FW_mic = [float(M[-3:]) for M in microturbulences[0]]
                VT = find_nearest(FW_mic,vmic)
                VT = 'VT000'[:-len(str(int(np.rint(VT))))]+str(int(np.rint(VT)))   
            except: 
                print('You should use a suitable microturbulence. Not found the correct VT')
                return None
        else: VT = vmic
        
        #######################Considered elemens abundances #######################
        elems = 'C','N','O','H_He'
        abundances = {}
        for elem,E in zip(elems,[E_C, E_N, E_O, E_H_He]):
            metallicities = pd.read_csv(dir().modeldir+'grid_components/abundances_'+elem+'.dat', header=None,delim_whitespace=True)[1] # gbat grid!
            #metallicities = pd.read_csv('../grid_components/abundances_'+elem+'_old.dat', header=None,delim_whitespace=True)[1] #OLD MODELS
            if E>0: E -= 12
            abundances[elem] = find_nearest(metallicities, E)
        print('Abundance Nitrogen: ',abundances['N']+12)
        ############ Finding the best model in a given list if not given ###########
        if 'CNO' not in model:
            model =  model+'CNO'
            
         ############## Finding vsin and vmac in a list if not given ##############
        
        ######################### Generates line spectrum #########################
        self.abundances, self.vmic = abundances,VT
        self.synthetic = spectrum_generation(model,elems,self.abundances,self.vmic,**kwargs)
        FL_conv = np.zeros(len(self.synthetic[1]))
        
        ####################### Reading broadening sections #######################
        broadening_sections = pd.read_csv('line_broadening_divisions.txt',skiprows=1, header=None,delim_whitespace=True)
        
        ##########################Broadening the spectra ##########################
        for i in broadening_sections.index:
            m = broadening_sections[0][i],broadening_sections[1][i]
            l0 = np.mean(m)
            M = (self.synthetic[0]>m[0]) & (self.synthetic[0]<m[1])
            WAVE, FLUX = self.synthetic[0][M],self.synthetic[1][M]
            #WV, FL = lb.final_resolution(WAVE, FLUX-1, l0, vsini, vmac, R, beta=beta)
            broadening_class = lb.Broadening_profile(WAVE,l0,vsini,beta,vmac,R)
            WV, FL = broadening_class.convolution(FLUX-1)
        
            f1 = interp1d(WV, FL, fill_value="extrapolate")
            start_wave,stop_wave = WAVE[0],WAVE[-1]
            X = np.arange(start=start_wave,stop=stop_wave+1e-3,step=np.diff(WAVE)[0])
            Y = f1(X)
            FL_conv[M] += Y
        FL_conv = FL_conv
        self.wave = self.synthetic[0]
        self.flux = FL_conv+1
        return None

############################# ADAPTED FROM spec.py ############################
def f_gaussian(x, sigma):
    return np.exp(-(x/sigma)**2/2)

def f_gaussian1(x, A, lam0, sigma):
    # A -> Amplitude;  lam0 -> center
    return A*np.exp(-(x - lam0)**2/(2*sigma**2)) + 1

def f_gaussian_general(x, *params):
    f = 0
    for i in range(len(params)//3):
        A,lam0,sigma = params[3*i],params[3*i+1],params[3*i+2]
        f += f_gaussian1(x, A, lam0, sigma)-1
    f += 1
    return f

def gaussian_convolved(x,A,lam0,sigma,vsini,beta,vmac,**kwargs):
    gaussian = f_gaussian1(x, A, lam0, sigma)
    broad_profile = lb.Broadening_profile(x,lam0,vsini,beta,vmac)
    final_profile = broad_profile.convolution(gaussian-1)
    return final_profile[1] + 1

def syn_spec_fitline(model,line_c,E_C=8.3,E_N=7.6,E_O=8.75,vmic=9,line=None, \
    width=3, tol=150., func='g', iter=3, fw3414=False, info=False, outfit=False,\
    plot=False,lim=None,p0=None, **kwargs):

    '''
    Function to fit a spectral line from a synthetic spectrum to a function.
    Different functions account for different lines depending on their 
    natural profile (e.g. metallic lines should be fitted with either a 
    gaussian, Voigt, rotational profiles). See fitline_readme.txt and 
    fitline_diagram.txt included with pyIACOB for more details.
    Adapted from function in spec.py

    Parameters
    ----------
    model : str or array-like
        - If str: expected to be the name of the model to compute. It has the form "T340g410He10Q130b10CNO" 
        
    line_c : float
        Sets the central wavelenght of the line to search and fit.
    
    E_C, E_N, E_O : float, optionals
        Abundances of the different elements. For He, it uses E_N. Defaults are 8.3,7.6,8.75 (reference solar values)
        If values are in format log[X/H]+12, it substracts 12
        It has the form log(X/H) + 12
        IF MODEL IS NOT STR IT IS NOT USED
        
    vmic : str or float/int, optional
        Microturbulence to be considered. It has to be in the pkl, and usually has the form 'VT007'. Default is 9 
        If type is not str, it looks for the closest value and write the corresponding format.        
    
    line : float, optional
        Line to be studied. If None, it takes the first line in line_c. Default is None

    width : int, optional
        Sets the width in [AA] where the line fits well in. Default is 3.

    tol : int, optional
        Sets the tolerance [km/s] to shifting the spectrum in order to fit the line.

    func : str, optional
        Choose the function to fit the line:
        'g' Gaussian (default); 'l' Lorentzian; 'v' Voigt; 'r' Rotational.
        'vr_H' Voigt profile with rotation for H lines.
        'vr_Z' Voigt profile with rotation for metallic lines.
        'vrg_H' Voigt profile with rotation and extra gaussian wings for H lines.
        'vrg_Z' Voigt profile with rotation and extra gaussian wings for metallic lines.

    iter : int, optional
        Number of iterations to optimize window width. Default is 3.

    fw3414 : boolean, optional
        If 'True', it will output the FW difference at 3/4 - 1/4 of the line depth.
        Default is 'False'.

    info : boolean, optional
        If 'True', it will print information for each line fitting.
        Default is 'False'.

    outfit : boolean, optional
        If 'True', it will output the fitting as an array.
        Default is 'False'.

    plot : boolean, optional
        If 'True', it will create and show the plots.
        
    lim : array-like, optional
        Limits to move the values to fit. Default is None    
    
    p0 : array-like, optional 
        Values to start to iterate. Default is None

    **kwargs : args for synthetic_model

    Returns
    -------
    Dictionary containing the results, fitting parameters and optionally the
    original line, normalized flux and fitted flux.

    Notes
    -----
    Emission lines are excluded, see "Filtering emission lines" section.
    Returns: Parameters from the fitted line and a last value containing data
    with input wavelenght and flux, plus the flux normalized, the flux of the
    fitted line, and the fitting parameters.
    '''
    if type(model) == str:
        print('str case')
        #######################Considered elemens abundances #######################
        #E_C, E_N, E_O, vmic = Model[star][4:]
        if E_C>0 or E_N>0 or E_O>0:
            E_C,E_N,E_O = E_C-12,E_N-12,E_O-12 
        E_H_He = E_N
        ##### Reading possible microturbulence and finding the closest to the given #####
        if type(vmic) != str:
            try:
                microturbulences = pd.read_csv(dir().modeldir+'grid_components/microturbulences.dat', header=None,delim_whitespace=True)      
                FW_mic = [float(M[-3:]) for M in microturbulences[0]]
                VT = find_nearest(FW_mic,vmic)
                VT = 'VT000'[:-len(str(int(np.rint(VT))))]+str(int(np.rint(VT)))   
            except: 
                print('You should use a suitable microturbulence. Not found the correct VT')
                return None
        else: VT = vmic
        #######################Considered elemens abundances #######################
        elems = 'C','N','O','H_He'
        #abundances = {}
        for elem in elems:
            metallicities = pd.read_csv(dir().modeldir+'grid_components/abundances_'+elem+'.dat', header=None,delim_whitespace=True)[1]
            #abundances[elem] = find_nearest(metallicities, locals()['E_'+elem]-12)
            
        #wave,flux = spectrum_generation(model,elems,abundances,VT)
        wave,flux = synthetic_model(None,E_C,E_N,E_O,vmic,model,**kwargs)
    else: 
        wave,flux = model
    
    fitsol = {'sol': 0, 'line': None, 'EW': None,'FWHM': None, 'depth': None, 'lines': None}

    #============================== Parameters =============================
    # Catch line input containing more than one line
    if line==None and type(line_c) == float:
        line = line_c
    elif line==None: 
        print('You have to add a central line')
        return fitsol
    if type(line_c) != float:
        line_c = np.array(line_c)
    #dlamb = line/self.resolution

    # Maximum shift between the minimum of the fitted line and the tabulated value
    tol_aa = float(tol)*(line)*1000/cte.c  # Changes km/s to angstroms

    # Maximum FWHM allowed (should be up to 18 for H lines)
    FWHM_max = 17

    #=========== Dictionary and boundary limits for the functions ==========
    fit_dic = {'g': ('A','lam0','sigma'),
               'l': ('A','lam0','gamma','y'),
               'v': ('A','lam0','sigma','gamma','y'),
               'r': ('A','lam0','sigma','vsini'),
              'vr': ('A','lam0','sigma','gamma','vsini','y'),
             'vrg': ('A','lam0','sigma','gamma','vsini','A2','sigma2','y'),
             'gc': ('A','lam0','sigma','vsini','beta','vmac')}

    # Fitting function: Gaussian | A,lam0,sig
    if func == 'g' and type(line_c) is float:
        fitfunc = f_gaussian1
        if lim==None:
            bounds  = ([-1,line-tol_aa,0],
                       [ 1,line+tol_aa,3]) # 3/6 for narrow/broad lines (default 6)
        else: bounds = lim
        
    elif func == 'g':
        
        fitfunc = f_gaussian_general
        
        if lim == None:  # We have to include deffault bounds. 
            bounds  = ([-1,line_c[0]-tol_aa,0,-1,line_c[1]-0.5,0],
                       [ 1,line_c[0]+tol_aa,3,0,line_c[1]+0.5,2]) # 3/6 for narrow/broad lines (default 6)
        else: bounds  = lim
    elif func == 'gc':
        fitfunc = gaussian_convolved #(x,A,lam0,sigma,vsini,beta,vmac,**kwargs)
        
        if lim == None:  # We have to include deffault bounds. 
            bounds  = ([-1,line_c-tol_aa,0,0,0,0],
                       [ 1,line_c+tol_aa,3,200,1.5,200]) # 3/6 for narrow/broad lines (default 6)            
        else: bounds  = lim    
    # Fitting function: Lorentzian | A,lam0,gamma,y



    #============================= Line fitting ============================
    i = 0; width_i = width
    p0_considered_def = []
    while i < iter:
        # Extracting the window of the spectrum
        if type(line_c) == float: window = (wave >= line-width_i/2) & (wave <= line+width_i/2)
        else:  window = (wave >= min(line_c)-width_i/2) & (wave <= max(line_c)+width_i/2)

        flux = flux[window]  # Selectingf the windows around the line
        wave = wave[window]
        dlam_mean = np.mean(np.diff(wave))

        # Final normalization of the iteration
        flx,wv = flux,wave
        flux_norm_i = flux
        #========================= Fitting the line ========================ç
        
        try:    
            if func=='g':
                presented_lines = []
                for k in range(len(p0)//3):
                    c,w = p0[3*k+1], 2.335*p0[3*k+2]
                    msk = (wave>(c-w/4)) & (wave<(c+w/4))
                    if abs(np.mean(flux[msk])-1)>2*np.std(flx):
                        presented_lines.append(True),presented_lines.append(True),presented_lines.append(True)
                    else:
                        print(c,abs(np.mean(flux[msk])-1),2*np.std(flx))
                        presented_lines.append(False),presented_lines.append(False),presented_lines.append(False)

                bounds_considered = np.array(bounds[0])[presented_lines], np.array(bounds[1])[presented_lines]
                p0_considered = np.array(p0)[presented_lines]
            else: 
                p0_considered = np.array(p0)
                bounds_considered = np.array(bounds[0]), np.array(bounds[1])

            popt_i = curve_fit(fitfunc, wave, flux_norm_i, bounds=bounds_considered, p0=p0_considered)[0]
            print('EW =',round(sum((flux_norm_i-1)*-10),2))
            p0_considered_def = p0_considered.copy()
            
            flux_fit_i = fitfunc(wave, *popt_i)
            # Calculate the empirical approximate FWHM
            medval = (max(flux_fit_i) + min(flux_fit_i))/2
            if popt_i[0]<0:
                medpos = [np.where(flux_fit_i <= medval)[0][value] for value in (0,-1)]
            else:
                medpos = [np.where(flux_fit_i >= medval)[0][value] for value in (0,-1)]
            FWHM = round(wave[medpos[1]] - wave[medpos[0]], 2)
            # Checking step results
            # The min FWHM will be defined by either 3 times the dlam, or 3/4 of the
            # minimum theoretical FWHM given the input line and resolution
            #FWHM_min = np.max([3*dlam_mean,3/4*dlamb])
            FWHM_min = 3*dlam_mean
            
            if FWHM_min < FWHM < FWHM_max:
                flux_norm = flux_norm_i; flux_fit = flux_fit_i; popt = popt_i; width = width_i
                i = i + 1; width_i = FWHM*7

                # We compute the gaussian corresponding to the studied line. 
                flux_fit_EW = f_gaussian1(wave, *popt_i[:3])
            elif FWHM < FWHM_min:
                print('WARNING: FWHM(%.1f) < minimum FWHM for %.3fA' % (FWHM,line))
                break
            elif FWHM > FWHM_max:
                print('WARNING: FWHM(%.1f) > maximum FWHM for %.3fA ' % (FWHM,line))
                break
            
        except:
            
            
            try: 
                popt_i = curve_fit(fitfunc, wave, flux_norm_i)
                print('Check exception try ',popt_i,fitfunc)
                flux_fit_i = fitfunc(wave, *popt_i)

                
                
                # Calculate the empirical approximate FWHM
                medval = (max(flux_fit_i) + min(flux_fit_i))/2
                if popt_i[0]<0:
                    medpos = [np.where(flux_fit_i <= medval)[0][value] for value in (0,-1)]
                else:
                    medpos = [np.where(flux_fit_i >= medval)[0][value] for value in (0,-1)]
                FWHM = round(wave[medpos[1]] - wave[medpos[0]], 2)
                # Checking step results
                # The min FWHM will be defined by either 3 times the dlam, or 3/4 of the
                # minimum theoretical FWHM given the input line and resolution
                #FWHM_min = np.max([3*dlam_mean,3/4*dlamb])
                FWHM_min = 3*dlam_mean
                if FWHM_min < FWHM < FWHM_max:
                    flux_norm = flux_norm_i; flux_fit = flux_fit_i; popt = popt_i; width = width_i
                    i = i + 1; width_i = FWHM*7

                    # We compute the gaussian corresponding to the studied line. 
                    flux_fit_EW = f_gaussian1(wave, *popt_i[:3])
                elif FWHM < FWHM_min:
                    print('WARNING: FWHM(%.1f) < minimum FWHM for %.3fA' % (FWHM,line))
                    break
                elif FWHM > FWHM_max:
                    print('WARNING: FWHM(%.1f) > maximum FWHM for %.3fA ' % (FWHM,line))
                    break
            except:             
                plt.figure()
                plt.plot(wave, flux_norm_i)
                print('Not feasable. Try including p0')
                break
    fitsol['lines'] = p0_considered_def
    print('Considered lines for ',line,fitsol['lines'])
    
    #======================== Checking final results =======================
    #window = (self.wave >= line-width/2.) & (self.wave <= line+width/2.)
    
    if type(line_c) == float: window2 = (wave >= line-width/2) & (wave <= line+width/2)
    else:  window2 = (wave >= min(line_c)-width/2) & (wave <= max(line_c)+width/2)

    flux = flux[window2]
    wave = wave[window2]

    if i == 0:
        if info is True:
            print('Problem in spectrum %s' % model)
            print('Line %sA could not be fitted or does not exist.\n' % line)

    if popt_i[0]<0:
        line_f = wave[flux_fit == min(flux_fit)][0]
    else:
        line_f = wave[flux_fit == max(flux_fit)][0]
    
    if abs(line - line_f) > tol_aa:
        if info is True:
            print('Line %sA found outside tolerance.\n' % line)
        return fitsol

    line_f = round(line_f, 3)

    #=========================== Calculate the EW ==========================
    # stackoverflow.com/questions/34075111/calculate-equivalent-width-using-python-code
    # When emission is considered abs should be removed and 1-flux_fit -> flux_fit-1
    # Compute the EW from the flux of the gaussian FROM the studied line. 
    EW = 0.5*abs(fsum((wave[wl-1] - wave[wl])*((1 - flux_fit_EW[wl-1]) +
         (1 - flux_fit_EW[wl])) for wl in range(1, len(flux_fit_EW))))
    EW = round(1000*EW,2)

    #====================== Calculate the final FWHM =======================
    medval = (max(flux_fit) + min(flux_fit))/2
    if popt_i[0]<0:
        medpos = [np.where(flux_fit <= medval)[0][value] for value in (0,-1)]
        try:
            l_val = np.interp(medval, [flux_fit[medpos[0]], flux_fit[medpos[0]-1]],
                [wave[medpos[0]], wave[medpos[0]-1]])
        except:
            l_val = wave[medpos[0]]
        try:
            r_val = np.interp(medval, [flux_fit[medpos[1]], flux_fit[medpos[1]+1]],
                [wave[medpos[1]], wave[medpos[1]+1]])
        except:
            r_val = wave[medpos[1]]

        FWHM = round(r_val - l_val, 3)
    else:
        medpos = [np.where(flux_fit >= medval)[0][value] for value in (0,-1)]
        try:
            l_val = np.interp(medval, [flux_fit[medpos[0]], flux_fit[medpos[0]-1]],
                [wave[medpos[0]], wave[medpos[0]-1]])
        except:
            l_val = wave[medpos[0]]
        try:
            r_val = np.interp(medval, [flux_fit[medpos[1]], flux_fit[medpos[1]+1]],
                [wave[medpos[1]], wave[medpos[1]+1]])
        except:
            r_val = wave[medpos[1]]

        FWHM = round(r_val - l_val, 3)            
    #======= Calculate width difference at 3/4-1/4 of the line depth =======
    if fw3414 is True:
        try:
            lowval = (max(flux_fit) + 3*min(flux_fit))/4
            uppval = (3*max(flux_fit) + min(flux_fit))/4

            medpos = []; FW34_14 = []
            for par,val in zip(['FW14_Hb','FW34_Hb'],[lowval,uppval]):
                medpos_i = [np.where(flux_fit <= val)[0][value] for value in (0,-1)]
                medpos.append(medpos_i)
                try: l_val = np.interp(val,[flux_fit[medpos_i[0]],flux_fit[medpos_i[0]-1]],
                                               [wave[medpos_i[0]],wave[medpos_i[0]-1]])
                except: l_val = wave[medpos_i[0]]
                try: r_val = np.interp(val,[flux_fit[medpos_i[1]],flux_fit[medpos_i[1]+1]],
                                              [wave[medpos_i[1]],wave[medpos_i[1]+1]])
                except: r_val = wave[medpos_i[1]]

                FW34_14.append(round(r_val-l_val, 3))

            FW34_14 = FW34_14[1]-FW34_14[0]

        except:
            print('Problem calculating the FW at 1/4 and 3/4 of the Hb line.')
            FW34_14 = np.nan
    else:
        FW34_14 = np.nan

    #======================= Calculate the line depth ======================
    depth = round(1 - min(flux_fit_EW), 3)

    #================================ Plot =================================
    if plot is True:
        fig, ax = plt.subplots()

        ax.tick_params(which='minor', length=2),ax.tick_params(which='major', length=5)   
        ax.xaxis.set_minor_locator(AutoMinorLocator(5))
        ax.yaxis.set_minor_locator(AutoMinorLocator(5))


        if fw3414 is True and FW34_14 != np.nan:
            print(medpos,lowval,uppval)
            ax.plot([wave[medpos[0][0]],wave[medpos[0][1]]],[lowval,lowval],'gray',linestyle='--',lw=.5)
            ax.plot([wave[medpos[1][0]],wave[medpos[1][1]]],[uppval,uppval],'gray',linestyle='--',lw=.5)

        ax.plot(wave, flux, c='orange', lw=.5)
        ax.plot(wave, flux_fit, c='g', lw=.5)
        ax.plot(wave, flux_norm, c='b', lw=.5)
        ax.set_title('EW: %d | FWHM: %.2f | FW34-14: %.2f' %
            ( EW,FWHM,FW34_14), fontsize=8)
        #ax.set_yticks([])
        ax.set_xlabel('$\lambda$ $[\AA]$', size=13)
        ax.set_ylabel('Normalized flux', size=13)
        ax.tick_params(direction='in', top='on', labelright=True, labelleft=False)
        ax.figure.subplots_adjust(top=.9, bottom=.12, right=.88, left=.08)
        ax.set_xlim(max(line-3*FWHM,plt.xlim()[0]),min(line+3*FWHM,plt.xlim()[1]))
        
    plt.show(block=False)

    #=======================================================================
    #================= Packing the results in a dictionary =================
    
    fitsol = {'sol':1, 'line':line_f,'EW':EW, 'FWHM':FWHM, 'FW34_14':FW34_14,
              'depth':depth,'lines': p0_considered_def}
    if info is True:
        print(line,popt_i)
    for f_par,par in zip(fit_dic[func.split('_')[0]], popt):
        fitsol[f_par] = round(par, 3)

    if outfit == True:
        fitsol['wave'] = wave
        #fitsol['flux'] = flux
        fitsol['flux_norm'] = flux_norm
        fitsol['flux_fit'] = flux_fit
    return fitsol

if __name__ == "__main__": 
    model_example = synthetic_model(None,E_C=8.26,E_N=7.65,E_O=8.8,vmic=12,model='T330g410He10Q1400b10CNO',\
                                folder_original_pkl=dir().fastwinddir+dir().fastwinddir.split('/')[-2]+'_pkl/original/')
