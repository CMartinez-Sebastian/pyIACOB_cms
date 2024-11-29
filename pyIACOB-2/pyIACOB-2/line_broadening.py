#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: C. Martinez-Sebastian

GOAL:
    Compute different components that broad the spectrum.    
"""

from __future__ import print_function, division
import numpy as np
import scipy.special as spc
from scipy import stats
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import tools as tls
import scipy.constants as cte


# In[rotational]
class rot:
    def __init__(self, vsini, beta):
        """
        Calculate the broadening profile due to rotation.
        Parameters
        ----------
        vsini : float
            Projected rotation speed of the star [km/s]
        beta : float
            Limb-darkening coefficient
        """
        self.vc = vsini / 299792.458
        self.bet = beta

    def rotational(self, dl, l0): # Correction to Carroll 1933
        """
        Calculates the broadening profile.
        Parameters
        ----------
        dl : array
            'Delta wavelength': The distance to the reference point in
            wavelength space [A].
        l0 : array
            The reference wavelength [A].
        Returns
        -------
        Broadening profile : array
            The broadening profile according to Gray. 
        """
        self.EPS = 1. / (1. + self.bet)
        self.dl0 = self.vc * l0
        self.c_b = 1. / (1. + 2. * self.bet / 3.)
        self.c_a = 1. / self.dl0
        self.c1 = 2. / np.sqrt(np.pi)
        self.c2 = self.bet / 2.
        self.x = dl/self.dl0
        self.mask = abs(self.x)<=1
        #print(self.mask)
        ##indi = np.where(np.abs(x) < 1.0)[0]
        ##result[indi] = self.c1 * \
        ##    np.sqrt(1. - x[indi]**2) + self.c2*(1. - x[indi]**2)
        #result = np.zeros(len(self.x))
        #result[self.mask] = self.c_b *  (self.c1 * np.sqrt(1. - self.x[self.mask]**2) + \
        #                      self.c2 * (1. - self.x[self.mask]**2)) * self.c_a
        result = self.c_b *  (self.c1 * np.sqrt(1. - self.x[self.mask]**2) + \
                             self.c2 * (1. - self.x[self.mask]**2)) * self.c_a
        
        # Correct the normalization for numeric accuracy
        # The integral of the function is normalized, however, especially in the case
        # of mild broadening (compared to the wavelength resolution), the discrete
        # broadening profile may no longer be normalized, which leads to a shift of
        # the output spectrum, if not accounted for.
        return result

def rotBroad(wvl, flux, vsini, l0, beta=1.5, edgeHandling="firstlast", ret="full"):
    """
    Apply rotational broadening to a spectrum.
    This function applies rotational broadening to a given
    spectrum using the formulae given in Gray's "The Observation
    and Analysis of Stellar Photospheres". It allows for
    limb darkening parameterized by the linear limb-darkening law.
    The `edgeHandling` parameter determines how the effects at
    the edges of the input spectrum are handled. If the default
    option, "firstlast", is used, the input spectrum is internally
    extended on both sides; on the blue edge of the spectrum, the
    first flux value is used and on the red edge, the last value
    is used to extend the flux array. The extension is neglected
    in the return array. If "None" is specified, no special care
    will be taken to handle edge effects.
    .. note:: Currently, the wavelength array as to be regularly
              spaced.
    Parameters
    ----------
    wvl : array
        The wavelength array [A]. Note that a
        regularly spaced array is required.
    flux : array
        The flux array.
    vsini : float
        Projected rotational velocity [km/s].
    l0 : float
        Central wavelength desired to compute the rotational shape. 
    beta : float, optional
        Limb-darkening coefficient (0-1.5). Default is 1.5
    edgeHandling : string, {"firstlast", "None"}
        The method used to handle edge effects.
    ret : string, {"full", "original"}
        To choose if the returned array has the length of the convolution or the
        same as the original. 
    Returns
    -------
    wave : array
        Wavelength base array of the covolved line
    flux : array
        Flux array after covolving with rotational broadening.
    """
    # Check whether wavelength array is evenly spaced
    l1_or, l2_or = wvl[0], wvl[-1]
    sp = wvl[1::] - wvl[0:-1]
    if abs(max(sp) - min(sp)) > 1e-6:
        print("Input wavelength array is not evenly spaced. Please, use evenly spaced input array.")
    if vsini <= 0.0:
        print("vsini must be positive.")
    if (beta < 0) or (beta > 1.5):
        print("Linear limb-darkening coefficient, beta, should be '0 < beta < 1.5'. Adapt beta.")

    # Wavelength binsize
    dwl = wvl[1] - wvl[0]

    if edgeHandling == "firstlast":
        # Number of bins additionally needed at the edges
        binnu = int(np.floor(((vsini / 299792.458) * max(wvl)) / dwl)) + 1
        # Defined 'valid' indices to be returned
        # Adapt flux array
        front = np.ones(binnu) * flux[0]
        end = np.ones(binnu) * flux[-1]
        flux = np.concatenate((front, flux, end))
        # Adapt wavelength array
        front = (wvl[0] - (np.arange(binnu) + 1) * dwl)[::-1]
        end = wvl[-1] + (np.arange(binnu) + 1) * dwl
        wvl = np.concatenate((front, wvl, end))
    else:
        print("Edge handling method '" + str(edgeHandling) + "' currently not \
              supported. Choose ones of the valid edge handling methods")
    dl = wvl - l0
    gdl = rot(vsini, beta)
    g = gdl.rotational(dl, l0)
    result = np.convolve(flux, g/sum(g))
    
    L1, L2 = wvl[0],wvl[-1]
    if ret == 'full':
        return np.linspace(L1+gdl.x[gdl.mask][0]*gdl.dl0,L2+gdl.x[gdl.mask][-1]*gdl.dl0,len(result)), result
    elif ret == 'original':
        wave = np.linspace(L1+gdl.x[gdl.mask][0]*gdl.dl0,L2+gdl.x[gdl.mask][-1]*gdl.dl0,len(result))
        eps = np.mean(sp)/2
        return wave[(wave >= l1_or - eps) & (wave <= l2_or + eps)], result[(wave >= l1_or - eps) & (wave <= l2_or + eps)]

# In[macroturbulence]
class mac:
    def __init__(self, vmac):
        """
        Calculate the broadening profile due to macroturbulence.
        Parameters
        ----------
        vmac : float
            Effective macroturbulence speed in the star [km/s]
        """
        self.vmacc = vmac / 299792.458

    def macroturbulence(self, dl, l0): # Correction to Carroll 1933
        """
        Calculates the broadening profile.
        Parameters
        ----------
        dl : array
            'Delta wavelength': The distance to the reference point in
            wavelength space [A].
        l0 : array
            The reference wavelength [A]
        Returns
        -------
        Broadenign : array
            Broadening profile derived from Gray 
        """
        self.dM = l0 * self.vmacc
        self.c1 = 2. / np.sqrt(np.pi)
        x = dl / self.dM
        self.x = x[x>=0]
        result = self.c1 * (- self.x * np.sqrt(np.pi) + np.exp(-self.x**2) + \
                            self.x * np.sqrt(np.pi) * spc.erf(self.x)) / self.dM
        # As there is a degenerancy between A_R, A_T, vmac_R and vmac_T in
        # eq. A.3 of Simon-Diaz's thesis, we fix three of the parameters, leaving
        # vmac (a projected velocity in km/s) as the only parameter.
        result = np.concatenate([np.flip(result),result])
        return result
def macBroad(wvl, flux, vmac, l0, edgeHandling="firstlast", ret="full"):
    """
    Apply macroturbulence broadening to a spectrum.
    This function applies macroturbulence broadening to a given
    spectrum using the formulae given in Gray's "The Observation
    and Analysis of Stellar Photospheres". It has been implemented 
    in the form of Simon-Diaz's thesis.
    The `edgeHandling` parameter determines how the effects at
    the edges of the input spectrum are handled. If the default
    option, "firstlast", is used, the input spectrum is internally
    extended on both sides; on the blue edge of the spectrum, the
    first flux value is used and on the red edge, the last value
    is used to extend the flux array. The extension is neglected
    in the return array. If "None" is specified, no special care
    will be taken to handle edge effects.
    .. note:: Currently, the wavelength array as to be regularly
              spaced.
    Parameters
    ----------
    wvl : array
        The wavelength array [A]. Note that a
        regularly spaced array is required.
    flux : array
        The flux array.
    vmac : float
        Effective macroturbulence speed in the star [km/s]
    l0 : float
        Central wavelength desired to compute the rotational shape
    edgeHandling : string, {"firstlast", "None"}
        The method used to handle edge effects.
    ret : string, {"full", "original"}
        To choose if the returned array has the length of the convolution or the
        same as the original. 
    Returns
    -------
    wave : array
        Wavelength base array of the covolved line
    flux : array
        Flux array after covolving with macroturbulence broadening.
    """
    # Check whether wavelength array is evenly spaced
    l1_or, l2_or = wvl[0], wvl[-1]
    sp = wvl[1::] - wvl[0:-1]
    if abs(max(sp) - min(sp)) > 1e-6:
        print("Input wavelength array is not evenly spaced. Please, use evenly spaced input array.")
    if vmac <= 0.0:
        print("vmac must be positive.")

    # Wavelength binsize
    dwl = wvl[1] - wvl[0]
    if edgeHandling == "firstlast":
        # Number of bins additionally needed at the edges
        binnu = int(np.floor(((vmac / 299792.458) * max(wvl)) / dwl)) + 1
        # Defined 'valid' indices to be returned
        # Adapt flux array
        front = np.ones(binnu) * flux[0]
        end = np.ones(binnu) * flux[-1]
        flux = np.concatenate((front, flux, end))
        # Adapt wavelength array
        front = (wvl[0] - (np.arange(binnu) + 1) * dwl)[::-1]
        end = wvl[-1] + (np.arange(binnu) + 1) * dwl
        wvl = np.concatenate((front, wvl, end))
    else:
        print("Edge handling method '" + str(edgeHandling) + "' currently not \
              supported. Choose ones of the valid edge handling methods")
    dl = wvl - l0
    mdl = mac(vmac)
    m = mdl.macroturbulence(dl, l0)
    result = np.convolve(flux, m/sum(m))#, 'same')
    
    L1, L2 = wvl[0],wvl[-1]
    if ret == 'full':
        return np.linspace(L1-mdl.x[-1]*mdl.dM,L2+mdl.x[-1]*mdl.dM,len(result)), result
    elif ret == 'original':
        wave = np.linspace(L1-mdl.x[-1]*mdl.dM,L2+mdl.x[-1]*mdl.dM,len(result))
        eps = np.mean(sp)/2
        return wave[(wave >= l1_or - eps) & (wave <= l2_or + eps)], result[(wave >= l1_or - eps) & (wave <= l2_or + eps)]


# In[instrumental]
class ins:
    def __init__(self, Res):
        """
        Calculate the broadening profile due to the instrument.
        Parameters
        ----------
        Res : float
            Instrumental resolution
        """
        self.R = Res
    def instrumental(self, l0, stp = 'Default'):
        """
        Calculate the broadening profile 
        ----------        
        l0 : array
            The reference wavelength [A]
        stp : float, optional
            Step desired to sample the gaussian of the instrument. If none is 
            given, it is computed as 1/10 of the FWHM.
        Returns
        -------
        Broadened : 
            Function to convolve with.

        """
        self.dl = l0 / self.R
        # We create an array to cover all the instrumental points, but with an
        # step size determined by the instrumental resolution. 
        if stp == 'Default':
            self.stp = self.dl/10
        else:
            self.stp = stp
        self.l = np.arange(-3.5*self.dl,3.5*self.dl,self.stp)
        self.sig = self.dl / (2 * np.sqrt(2 * np.log(2)))
        result = 1./(np.sqrt(2.*np.pi)*self.sig)*np.exp(-np.power((self.l)/self.sig, 2.)/2)
        return result
def insBroad(wvl, flux, R, l0, edgeHandling="firstlast", ret="full",**kwargs):
    """
    Apply instrumental broadening to a spectrum.
    This function applies instrumental broadening to a given
    spectrum using a gaussian profile. 
    The `edgeHandling` parameter determines how the effects at
    the edges of the input spectrum are handled. If the default
    option, "firstlast", is used, the input spectrum is internally
    extended on both sides; on the blue edge of the spectrum, the
    first flux value is used and on the red edge, the last value
    is used to extend the flux array. The extension is neglected
    in the return array. If "None" is specified, no special care
    will be taken to handle edge effects.
    .. note:: Currently, the wavelength array as to be regularly
              spaced.
    Parameters
    ----------
    wvl : array
        The wavelength array [A]. Note that a regularly spaced array is required.
        In this case, the step size should coincide with this of the instrumental 
        function. If it is not the case, it is interpolated to the same resolution
    flux : array
        The flux array.
    R: float
        Instrumental resolution
    l0 : float
        Central wavelength desired to compute the rotational shape
    edgeHandling : string, {"firstlast", "None"}
        The method used to handle edge effects.
    ret : string, {"full", "original"}
        To choose if the returned array has the length of the convolution or the
        same as the original. 
    Returns
    -------
    wave : array
        Wavelength base array of the covolved line
    flux : array
        Flux array after covolving with instrumental broadening.
    """
    # Check whether wavelength array is evenly spaced
    l1_or, l2_or = wvl[0], wvl[-1]
    sp = wvl[1::] - wvl[0:-1]
    if abs(max(sp) - min(sp)) > 1e-6:
        print("Input wavelength array is not evenly spaced. Please, use evenly spaced input array.")
    if R <= 0.0:
        print("vmac must be positive.")

    # Wavelength binsize
    dwl = wvl[1] - wvl[0]
    if edgeHandling == "firstlast":
        # Number of bins additionally needed at the edges
        # binnu = int(np.floor(((vmac / 299792.458) * max(wvl)) / dwl)) + 1     # Hay que cambiar esyto!!!
        binnu = int(np.floor(max(wvl) / dwl / R)) + 1
        # binnu = 150
        # Defined 'valid' indices to be returned
        # Adapt flux array
        front = np.ones(binnu) * flux[0]
        end = np.ones(binnu) * flux[-1]
        flux = np.concatenate((front, flux, end))
        # Adapt wavelength array
        front = (wvl[0] - (np.arange(binnu) + 1) * dwl)[::-1]
        end = wvl[-1] + (np.arange(binnu) + 1) * dwl
        wvl = np.concatenate((front, wvl, end))
    else:
        print("Edge handling method '" + str(edgeHandling) + "' currently not \
              supported. Choose ones of the valid edge handling methods")
    #dl = wvl - l0
    idl = ins(R)
    i = idl.instrumental(l0,**kwargs)
    if abs(np.mean(np.diff(wvl))-idl.stp)>0.1*idl.stp:
        f1 = interp1d(wvl,flux, kind='linear')  
        wvl = np.arange(start=wvl[0],stop=wvl[-1],step=idl.stp)
        flux = f1(wvl) 
    result = np.convolve(flux, i/sum(i))
    L1, L2 = wvl[0],wvl[-1]
    if ret == 'full':
        return np.linspace(L1+idl.l[0],L2+idl.l[-1],len(result)), result
    elif ret == 'original':
        wave = np.linspace(L1+idl.l[0],L2+idl.l[-1],len(result))
        eps = np.mean(sp)/2
        return wave[(wave >= l1_or - eps) & (wave <= l2_or + eps)], result[(wave >= l1_or - eps) & (wave <= l2_or + eps)]

# In[Total_convolution]

def final_resolution(wave, flux, l0, vsini, vmac, R, beta=1.5, stp ='Default', **kwargs):
    """
    Apply macroturbulence, rotational and instrumental broadening to a given spectrum,
    and interpolates to the given resolution (marked with parameter stp). 
    
    Parameters
    ----------
    wave : array,
        Original wavelength.
    flux : array,
        Original flux (before broadening.
    l0 : float,
        Central wavelength to consider the functions to broad.
    vsini : float,
        Projected rotational velocity (in km/s).
    vmac : float,
        Macroturbulence velocity (in km/s).
    R : float,
        Instrumental resoltuon used.
    beta : float, optional
        DESCRIPTION. The default is 1.5.
    stp : float, optional
        DESCRIPTION. The default is 'Default'.
    **kwargs : Non-specified common function arguments. 
        edgeHandling="firstlast", ret="full"
        
    Warnings
    --------
        As the broadening functions are dependent of the central wavelength 
        considered (which chnnges depending on the line), it is advisable to 
        divide a complete spectrum in smaller pieces to increase the accuracy of 
        the central wavelength 

    Returns
    -------
    wvl : array,
        Corresponding wavelength for the covolved flux.
    flux : array,
        Flux after convolutiuon with the same step size than at the begginng.

    """
    if stp=='Default':
        stp = l0 / (R * 10)
    f1 = interp1d(wave,flux, kind='linear')  
    wv = np.arange(start=wave[0],stop=wave[-1],step=stp)
    flx = f1(wv) 
    wvl, flux = rotBroad(wv, flx, vsini, l0, beta=beta, **kwargs)
    wvl, flux = macBroad(wvl, flux, vmac, l0, **kwargs)
    wvl, flux = insBroad(wvl, flux, R, l0, stp = stp, **kwargs)
    return wvl, flux
# In[Broadening_class]

class Broadening_profile:
    def __init__(self,wave,l0,vsini,beta,vmac,R=85000,stp=None,info=False):
        if stp==None:
            self.stp = np.mean(np.diff(wave))
        else: self.stp = stp
        self.wave = wave
        self.wv = np.arange(start=self.wave[0],stop=self.wave[-1]+self.stp/2,step=self.stp)
        if info: print('Step is ',self.stp)
        # ================= Creating broadening function ================
        # Computing the different broadening profiles
        self.rotational_profile = 1
        self.macroturbulence_profile = 1
        if vsini>1e-2:
            self.rot_profile = rot(vsini,beta)
            rotation_profile = self.rot_profile.rotational(self.wv-l0,l0)
            self.rotational_profile = rotation_profile/sum(rotation_profile)
            self.wvl_rot = np.linspace(-(len(self.rotational_profile)-1)/2*np.mean(np.diff(self.wv)),\
                                  (len(self.rotational_profile)-1)/2*np.mean(np.diff(self.wv)),len(self.rotational_profile))        
        if vmac>1e-2:    
            self.mac_profile = mac(vmac)
            macroturbulence_profile = self.mac_profile.macroturbulence(self.wv-l0,l0)
            
            
            self.macroturbulence_profile = macroturbulence_profile/sum(macroturbulence_profile)
            self.wvl_mac = np.linspace(-(len(self.macroturbulence_profile)-1)/2*np.mean(np.diff(self.wv)),\
                                  (len(self.macroturbulence_profile)-1)/2*np.mean(np.diff(self.wv)),len(self.macroturbulence_profile))
        self.ins_profile = ins(R)
        instrumental_profile = self.ins_profile.instrumental(l0,stp=self.stp)
        self.instrumental_profile = instrumental_profile/sum(instrumental_profile)
        self.wvl_ins = np.linspace(-len(self.instrumental_profile[:instrumental_profile.argmax()])*np.mean(np.diff(self.wv)),\
                              len(self.instrumental_profile[instrumental_profile.argmax():-1])*np.mean(np.diff(self.wv)),len(self.instrumental_profile))
        
        # Trying to covolve the profiles
        self.convolve_rot_mac = np.convolve(self.rotational_profile,self.macroturbulence_profile)
        self.convolve_rot_mac = self.convolve_rot_mac[self.convolve_rot_mac>0]
        self.wvl_rot_mac =  np.linspace(-len(self.convolve_rot_mac[:self.convolve_rot_mac.argmax()])*np.mean(np.diff(self.wv)),\
                              len(self.convolve_rot_mac[self.convolve_rot_mac.argmax():-1])*np.mean(np.diff(self.wv)),len(self.convolve_rot_mac))
        
        self.convolve_rot_mac_ins = np.convolve(self.convolve_rot_mac,self.instrumental_profile)
        self.wvl_rot_mac_ins =  np.linspace(-len(self.convolve_rot_mac_ins[:self.convolve_rot_mac_ins.argmax()])*np.mean(np.diff(self.wv)),\
                              len(self.convolve_rot_mac_ins[self.convolve_rot_mac_ins.argmax():-1])*np.mean(np.diff(self.wv)),len(self.convolve_rot_mac_ins))
        
    def convolution(self,flux):
        f1 = interp1d(self.wave,flux, kind='linear',fill_value="extrapolate")  
        self.flx = f1(self.wv)
        total_convolve = np.convolve(self.convolve_rot_mac_ins,self.flx)
        
        self.total_convolve = total_convolve
        self.total_convolve_cut = total_convolve[len(self.convolve_rot_mac_ins)//2-1:-len(self.convolve_rot_mac_ins)//2+1]
        self.total_wave = np.linspace((self.wv[0]+self.wv[1])/2-np.mean(np.diff(self.wv)),\
                                      (self.wv[-1]+self.wv[-2])/2+np.mean(np.diff(self.wv)),\
                                      len(self.total_convolve_cut))
        f1 = interp1d(self.total_wave,self.total_convolve_cut, kind='linear',fill_value="extrapolate")  
        self.flx_interp = f1(self.wave)         
        return self.wave,self.flx_interp 
    
    
    

# In[Total_convolution]

if __name__ == "__main__":
    import sys
    import matplotlib.pyplot as plt
    import FASTWIND_plotting as fw

    sys.path.append('/home/cms/Desktop/PhD/PROJECT/TESTS')
    from spectrum_generation_function import *

    #print(colored('RUNNING TEST', 'magenta', 'on_white',['bold']))
    print('\033[4mRUNNING TEST\033[0m')
    # Creating the synthetic spectra over which the test should be run
    elems = 'C','N','O','H_He'
    elems_abundance = {'C': -3.7, 'N': -4.2, 'O': -3.25, 'H_He': -4.2}
    wave, flux = spectrum_generation('T310g350He10Q1350b10CNO',elems,elems_abundance,'VT009')

    l0,vsini,vmac,R,beta,stp =5875,50,50,85000,1.5,0.005
    l0,vsini,vmac,R,beta,stp =5875,350,50,85000,1.5,0.005
    l0,vsini,vmac,R,beta,stp =5875,350,50,1000,1.5,0.005
    l0,vsini,vmac,R,beta,stp =5875,50,50,1000,1.5,0.005
    # ============================= Interpolating =============================
    if stp=='Default':
        stp = l0 / (R * 10)
    f1 = interp1d(wave,flux, kind='linear')  
    wv = np.arange(start=5800,stop=5950,step=stp)
    #wv = np.arange(start=6640,stop=6720,step=stp)
    flx = f1(wv) 
    # ============================= Default order =============================    
    wvl1,flux1 = final_resolution(wv, flx-1, l0, vsini, vmac, R, stp = stp)

    # ================= Trying to create a broadening function ================
    fig = plt.figure()
    fig.set_figheight(7)
    fig.set_figwidth(10)
    axs = {}
    fig.suptitle(r'$\lambda_0='+str(l0)+r';\,v\sin(i)='+str(vsini)+r';\,v_{\rm mac}='+str(vmac)+r';\,R='+str(R)+r'$')
    axs[0] = fig.add_subplot(3,2,1); axs[1] = fig.add_subplot(3,2,2)
    axs[2] = fig.add_subplot(3,2,3); axs[3] = fig.add_subplot(3,2,4)
    axs[4] = fig.add_subplot(3,1,3)
    fw.plot_style(list(axs.values()),size=12)
    
    axs[1].sharex(axs[0])
    axs[2].sharex(axs[0])
    axs[3].sharex(axs[0])
    axs[0].set_xlim(-10,10)
    # Computing the different broadening profiles
    
    ej = Broadening_profile(wv,l0,vsini,beta,vmac,R)

    
    rotation_profile = ej.rotational_profile
    wvl_rot= np.linspace(-(len(rotation_profile)-1)/2*np.mean(np.diff(wv)),\
                          (len(rotation_profile)-1)/2*np.mean(np.diff(wv)),len(rotation_profile))
    macroturbulence_profile = ej.macroturbulence_profile
    wvl_mac = np.linspace(-(len(macroturbulence_profile)-1)/2*np.mean(np.diff(wv)),\
                          (len(macroturbulence_profile)-1)/2*np.mean(np.diff(wv)),len(macroturbulence_profile))
    instrumental_profile = ej.instrumental_profile
    wvl_ins = np.linspace(-len(instrumental_profile)/2*np.mean(np.diff(wv)),\
                          len(instrumental_profile)/2*np.mean(np.diff(wv)),len(instrumental_profile))
    
    axs[0].plot(wvl_rot,rotation_profile);axs[0].set_title('Rotational profile',x=0.3,y=0.6)#;axs[0].set_xlim(-(vsini+5)/3e5*l0,(vsini+5)/3e5*l0)
    axs[1].plot(wvl_mac,macroturbulence_profile);axs[1].set_title('Macroturbulence profile',x=0.4,y=0.6)#;axs[1].set_xlim(-(vmac+5)/3e5*l0,(vmac+5)/3e5*l0)
    axs[2].plot(wvl_ins,instrumental_profile);axs[2].set_title('Instrumental profile',x=0.3,y=0.6)#;axs[2].set_xlim(-2*l0/R,2*l0/R)   
    
    wv_ej,flx_ej = ej.convolution(flx-1)
    
    wvl_broad = np.linspace(-len(ej.convolve_rot_mac_ins)/2*np.mean(np.diff(wv)),\
                          len(ej.convolve_rot_mac_ins)/2*np.mean(np.diff(wv)),len(ej.convolve_rot_mac_ins))
    axs[3].plot(wvl_broad,ej.convolve_rot_mac_ins);axs[3].set_title('Total convolution profile',x=0.4,y=0.6)
        
    axs[4].plot(wave,flux,label='Original FASTWIND spectrum')
    axs[4].plot(wvl1,flux1+1,label='Convolutioning EVERYTHING')
    axs[4].plot(wv_ej,flx_ej+1,'--',label='With broadening function')
    #axs[4].set_xlim(wv[0],wv[-1]),axs[4].legend()
    axs[4].set_xlim(5860,5890),axs[4].legend()
    
    print('Different EW:')
    mask_orig = (wave>6640) & (wave<6720) 
    mask_conv = (wvl1>6640) & (wvl1<6720) 
    mask_broad = (wv_ej>6640) & (wv_ej<6720) 
    print('EW_orig=',sum((flux[mask_orig]-1)*np.mean(np.diff(wave[mask_orig])))*-1e3)
    print('EW_conv=',sum(flux1[mask_conv]*np.mean(np.diff(wvl1[mask_conv])))*-1e3)
    print('EW_broad=',sum(flx_ej[mask_broad]*np.mean(np.diff(wv_ej[mask_broad])))*-1e3)
    fig.savefig('/home/cms/Desktop/PhD/HELPS/Mar/'+'lam'+str(l0)+'_vsini'+str(vsini)+'_vmac'+str(vmac)+'_R'+str(R)+'.pdf')