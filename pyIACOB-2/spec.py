from db import *

# Core packages
from math import fsum
import scipy.constants as cte
from scipy.special import wofz,erf
from scipy.optimize import curve_fit
from scipy.signal import convolve
from scipy.interpolate import interp1d
# Astro-packages
from astropy.time import Time
from astropy.stats import sigma_clip

# Plot packages
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator,LinearLocator
from line_broadening import Broadening_profile,final_resolution

import sys
sys.path.append('/home/cms/Desktop/PhD/PROJECT/TESTS')
from spectrum_generation_function import spectrum_generation
import formulas as fm

#plt.rc('xtick.minor', visible=True)
plt.rc('ytick', direction='in', right='on')
#plt.rc('ytick.minor', visible=True)

print('Reading spec...')
class spec():
    def __init__(self, spectrum, SNR=None, rv0=0, offset=0, txt=False, cut_edges=False):

        '''
        Parameters
        ----------

        spectrum : str
            Enter the input spectrum, either name(s) of the star(s), the fits files
            separated by coma, a .txt/.lst file containing the filenames, or '*'
            if you want to select all the fits files inside the working folder.

        SNR : str, optional
            If 'best' as input, it finds the best SNR spectrum for the given name.
            If 'bestHF' same as 'best' but prioritizing spectra from HERMES/FEROS.

        rv0 : float, optional
            Enter the radial velocity correction to apply to the spectrum in km/s.

        offset : float, optional
            Enter the offset in wavelength [A] of the spectrum to plot. Default is 0.

        txt : boolean, optional
            If True, it assumes spectrum from a two-columns file with wavelenght and flux
            with no header in the file. Default is False.

        cut_edges : boolean, optional
            If True, it cuts the edges of the spectrum where the flux is asymptotic.
            Default is False.
        '''

        if type(spectrum) == list:
            if len(spectrum) > 1:
                print('Error in spec(): More than one spectrum selected.\nExitting...')
                return None
            else: 
            	spectrum = spectrum[0]

        self.fullpath = findstar(spectrum, SNR=SNR)
        
        if self.fullpath is None:
            try: self.fullpath = [spectrum]
            except:
                if self.fullpath is None:
                    print('Error in spec(): No spectrum found.\nExitting...')
                    return None
            
        if len(self.fullpath) > 1:
            print(self.fullpath,len(self.fullpath))
            print('Error in spec(): More than one spectrum selected.\nExitting...')
            return None
        self.fullpath = self.fullpath[0]
        self.filename = self.fullpath.split(os.sep)[-1]
        self.id_star = self.fullpath.split(os.sep)[-1].split('_')[0]
        print(self.fullpath)
        self.resolution = int(re.split('(\d*\d+)',self.fullpath)[-2])

        self.offset = offset # Note, run self.waveflux to apply offset.

        self.rv0 = rv0 # Note, run self.waveflux to apply the correction.
        
        if txt == False:
            self.waveflux(cut_edges=cut_edges)
        elif txt == True:
            self.txtwaveflux(cut_edges=cut_edges)
        self.txt = txt
    def get_spc(self):

        '''
        Function to retrieve the spectral classification of the star from Simbad
        and add it to the class.
        '''

        try:
            query = Simbad.query_object(self.id_star)
            if query is None and 'HD' in self.id_star:
                new_id_star = self.id_star.replace('HD', 'HD ')
                query = Simbad.query_object(new_id_star)
            self.SpC = query['SP_TYPE'][0]
            #self.otypes = query['OTYPES'][0]

        except:
            print('Spectral classification could not be queried for %s' % self.id_star)
            self.SpC = ''


    def waveflux(self, lwl=None, rwl=None, width=0, helcorr='hel', cut_edges=False):

        '''
        Function to load or update the wavelenght and flux vectors and optionally apply
        an offset or a radial velocity correction if they are different from 0 in the
        class. It also adds the HJD, dlam to the class.

        Parameters
        ----------
        lwl : float, optional
            Sets the start wavelenght of the spectrum.

        rwl : float, optional
            Sets the end wavelenght of the spectrum.

        width : int, optional
            Sets the width in [AA] where the line fits well in. Default is 10.

        helcorr : str, optional
            If 'hel' as input (default), it applies the heliocentric correction.

        cut_edges : boolean, optional
            If True, it cuts the edges of the spectrum where the flux is asymptotic.
            Default is False.

        Returns
        -------
        In addition to update the class with new data, it returns the wavelength and flux
        vectors together with the HJD.
        '''

        # Retrieve the key values fron the fits header
        hdu = fits.open(self.fullpath)  # Open the fits image file
        hdu.verify('fix')               # Fix possible issues with the keywords
        header0 = hdu[0].header         # Read header of primary extension
        #print(header0)

        instrum = header0['INSTRUME']   # Instrument
        texp = header0['I-TEXP']          # Exposure time 
        lam0 = header0['CRVAL1']          # Get the wavelenght of the first pixel
        dlam = header0['CDELT1']          # Step of increase in wavelength
        pix0 = header0['CRPIX1']        # Reference pixel (generally 1, FEROS -49)
        spec_length = header0['NAXIS1'] # Length of the spectrum
        # Alternatively use len(hdu[0].data[0]) (NOT/MERCATOR) or len(hdu[0].data)

        # Correct Mercator CRVAL1 20101018-19:
        if any(bad in self.fullpath for bad in ['_20101018_','_20101019_']) and lam0 == 3763.9375:
            lam0 = 3763.61

        try:
            vbar = header0['I-VBAR'] # [km/s] Barycent. rv correction at midpoint
        #    vbar = header0['BVCOR']  # [km/s] Barycent. rv correction at midpoint | MERCATOR
        #    vbar = header0['VHELIO'] # [km/s] Barycent. rv correction at midpoint | NOT
        except:
            print('No helio/bary-centric correction applied to' + self.fullpath); vbar = 0

        self.vbar = vbar
        self.texp = texp

        try:
            hjd = header0['I-HJD']   # Heliocentric Julian date at midpoint
        except:
            hjd = Time(header0['DATE'], scale='utc').jd

        self.hjd = hjd

        try:
            self.snr = header0['I-SNR']  # SNR from header
        except:
            self.snr = np.nan

        try:
            self.SpC = header0['I-SPC']  # SpC from header
        except:
            self.SpC = ''

        # Make lists with wavelenght and flux for each spectrum
        if width >= 200:
            width = 200
            print('\nWARNING: Width value %f is too large, setting it to 200. ' %width)

        wave = lam0 + dlam*(np.arange(spec_length) - pix0 + 1)
        if '_log' in self.fullpath:
            wave = np.exp(wave)
        elif helcorr == 'hel' and not instrum == 'FEROS':
            wave = wave*(1 + 1000*vbar/cte.c)
        # Those with log and those from FEROS are already corrected from helcorr

        wave = wave*(1 - 1000*self.rv0/cte.c)

        wave = wave - self.offset

        try:
            flux = hdu[0].data[0]
        except:
            flux = hdu[0].data

        if cut_edges == True:
            if flux[1] < flux[0]:
                l_cut = next(pix for pix in range(len(wave)) if flux[pix+1] > flux[pix])
            else:
                l_cut = next(pix for pix in range(len(wave)) if flux[pix+1] < flux[pix])
            if flux[-2] < flux[-1]:
                r_cut = next(pix for pix in reversed(range(len(wave))) if flux[pix-1] > flux[pix])
            else:
                r_cut = next(pix for pix in reversed(range(len(wave))) if flux[pix-1] < flux[pix])

            wave = wave[l_cut:r_cut]
            flux = flux[l_cut:r_cut]

        if lwl != None and rwl != None:
            if wave[0] > lwl+dlam or wave[-1] < rwl-dlam:
                print('WARNING: Wavelenght limits outside spectrum wavelenght range.')
            flux = flux[(wave >= lwl-width/2.) & (wave <= rwl+width/2.)]
            wave = wave[(wave >= lwl-width/2.) & (wave <= rwl+width/2.)]

        if '_log' in self.fullpath:
            self.dlam = (wave[-1]-wave[0])/(len(wave)-1)
        else:
            self.dlam = dlam

        hdu.close()

        self.wave_0 = wave
        self.flux_0 = flux

        self.wave = wave
        self.flux = flux

        return wave, flux, hjd


    def txtwaveflux(self, lwl=None, rwl=None, width=0, cut_edges=False):

        '''
        Equivalent to spec.waveflux() but for spectra coming from ascii files.

        Note: The ascii files should be placed in a /ASCII/ subfolder inside data folder.

        Parameters
        ----------
        See help for spec.waveflux
        '''

        if self.filename.endswith('.fits'):
            print('Error in spec(): This is a fits file, use a .txt/.dat/.ascii or similar file.')
            return None
        data = findtable(self.filename, path=datadir+'ASCII/', format='basic')
        if data is None:
            data = findtable(self.filename,path=self.fullpath.replace(self.filename,''), format='basic')
        #print(data)
        if data.colnames[0].lower() in ['wave','wavelength','lambda','lamb','ang','angstroms']:
            wave = data[data.colnames[0]]
        else:
            print('Error in spec(): No wavelength column found in the firs column of the ascii file.')
            return None
        if data.colnames[1].lower() in ['flux','fluxes','norm_flux','flux_norm']:
            flux = data[data.colnames[1]]
        else:
            print('Error in spec(): No flux column found in the second column of the ascii file.')
            return None

        wave = wave*(1 - 1000*self.rv0/cte.c)

        wave = wave - self.offset

        if cut_edges == True:
            l_cut = next(pix for pix in range(len(wave)) if flux[pix+1] > flux[pix])
            r_cut = next(pix for pix in reversed(range(len(wave))) if flux[pix-1] > flux[pix])

            wave = wave[l_cut:r_cut]
            flux = flux[l_cut:r_cut]

        if lwl != None and rwl != None:
            dlam = (wave[-1]-wave[0])/(len(wave)-1)
            if wave[0] > lwl+dlam or wave[-1] < rwl-dlam:
                print('WARNING: Wavelenght limits outside spectrum wavelenght range.')
            flux = flux[(wave >= lwl-width/2.) & (wave <= rwl+width/2.)]
            wave = wave[(wave >= lwl-width/2.) & (wave <= rwl+width/2.)]

        self.dlam = (wave[-1]-wave[0])/len(wave)
        self.vbar = 0
        self.hjd = 0

        self.get_spc()

        self.wave_0 = wave
        self.flux_0 = flux

        self.wave = wave
        self.flux = flux
        return wave, flux, 0


    def fitline(self, line_c, line=None, width=15, tol=150., func='g', iter=3, fw3414=False, info=False,
        outfit=False, plot=False,lim=None,p0=None,normalization=True,sigma_tres=2.0,EW_err_k=False,iter_err=500,
        fit_range=None):

        '''
        Function to fit a spectral line to a function. Different functions account for
        different lines depending on their natural profile (e.g. metallic lines should be
        fitted with either a gaussian, Voigt, rotational profiles). See fitline_readme.txt
        and fitline_diagram.txt included with pyIACOB for more details.

        Parameters
        ----------
        line_c : float
            Sets the central wavelenght of the line to search and fit.

        width : int, optional
            Sets the width in [AA] where the line fits well in. Default is 15.

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

        fit_range : list
            Limits for the line to be fitted. If given, it give wave+flux in that range.
            Default is None
        Returns
        -------
        Dictionary containing the results, fitting parameters and optionally the
        original line, normalized flux and fitted flux.
        The EW in the dictionary is this of the first line in line_c (the one written in 'line_EW_studied')
        Notes
        -----
        Emission lines are excluded, see "Filtering emission lines" section.
        Returns: Parameters from the fitted line and a last value containing data
        with input wavelenght and flux, plus the flux normalized, the flux of the
        fitted line, and the fitting parameters.
        '''

        fitsol = {'sol': 0, 'line': None, 'RV_A': None, 'RV_kms': None, 'EW': None,
                  'EW_err': None, 'FWHM': None, 'depth': None, 'q_fit': None, 
                  'snr': None, 'lines': None, 'line_EW_studied' : None}

        #============================== Parameters =============================
        # Catch line input containing more than one line
        '''
        if type(line) == str and (',' in line or ' ' in line.strip()):
            print('Error in spec(): More than one line selected.\nExitting...')
            return fitsol
        '''
        if line==None and type(line_c) == float:
            line = line_c
        elif line==None and type(line_c) != float:
            line_c = np.array(line_c)
            line = line_c[0]
        elif line==None: 
            print('You have to add a central line')
            return fitsol
        if type(line_c) != float:
            line_c = np.array(line_c)
            if line==None: line = line_c[0]
        dlamb = line/self.resolution

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
        if func == 'gc' and type(line_c) is float:
            fitfunc = gaussian_convolved_general #(x,A,lam0,sigma,vsini,beta,vmac,**kwargs) 
            if lim == None:  # We have to include deffault bounds. # 3/6 for narrow/broad lines (default 6)  
                bounds  = ([-1,line_c-tol_aa,0,0,0,0],
                           [ 1,line_c+tol_aa,3,200,1.5,200]) # 3/6 for narrow/broad lines (default 6)             
            else: bounds  = lim
        elif func == 'gc':
            fitfunc = gaussian_convolved_general
            if lim == None:  # We have to include deffault bounds. 
                bounds  = ([*(-1*np.ones(len(line_c))),*(np.array(line_c-tol_aa)),0,0,0,0],
                           [*(1*np.ones(len(line_c))),*(np.array(line_c+tol_aa)),3,200,1.5,200])          
            else: bounds  = lim            
        
        if func == 'g' and type(line_c) is float:
            fitfunc = f_gaussian1
            if lim==None:
                bounds  = ([-1,line-tol_aa,0],
                           [ 0,line+tol_aa,3]) # 3/6 for narrow/broad lines (default 6)
            else: bounds = lim  
            
        elif func == 'g':
            
            fitfunc = f_gaussian_general
            arr = np.ones(len(line_c))
            if lim == None:  # We have to include deffault bounds. 
                #bounds  = ([-1,line_c[0]-tol_aa,0,-1,line_c[1]-0.5,0],
                #           [ 0.2,line_c[0]+tol_aa,3,0,line_c[1]+0.5,2]) # 3/6 for narrow/broad lines (default 6)
                
                bounds  = ([-1,-1,line_c[0]-tol_aa,line_c[1]-tol_aa,0,0],
                           [0.5,0.5,line_c[0]+tol_aa,line_c[1]+tol_aa,3,3]) # 3/6 for narrow/broad lines (default 6)
            else: bounds  = lim
            
        # Fitting function: Lorentzian | A,lam0,gamma,y
        elif func == 'l':
            fitfunc = f_lorentzian
            bounds  = ([-1,line-tol_aa, 0,1. ],
                       [ 0,line+tol_aa,10,1.01])

        # Fitting function: Voigt profile | A,lam0,sigma,gamma,y
        elif func == 'v':
            fitfunc = f_voigt
            bounds  = ([-10,line-tol_aa,0. ,0. ,1.  ], #'A' ~15 for As
                       [  0,line+tol_aa,7.5,8.5,1.01])

        # Fitting function: Rotational profile | A,lam0,sigma,vsini
        elif func == 'r':
            fitfunc = f_rot
            bounds  = ([.0,line-tol_aa,0. ,  1],
                       [.3,line+tol_aa,2.5,410])

        # Fitting function: Voigt x Rotational profile | A,lam0,sigma,gamma,vsini,y
        elif func == 'vr_H':
            fitfunc = f_voigtrot
            bounds  = ([-.3,line-tol_aa,0. ,0,  1,.0 ],
                       [ .0,line+tol_aa,4.5,8,410,.01])
        elif func == 'vr_Z':
            fitfunc = f_voigtrot
            bounds  = ([-.3,line-tol_aa,0. ,0,  1,.0 ],
                       [ .0,line+tol_aa,1.5,1,410,.01])

        # Fitting function: Voigt x Rotational + Gaussian profile
        # A,lam0,sigma,gamma,vsini,A2,sigma2,y
        elif func == 'vrg_H':
            fitfunc = f_vrg
            bounds  = ([-.5,line-tol_aa, 0, 0,  1,-.07,0,.0 ],
                       [ .0,line+tol_aa,10,10,410, .0 ,4,.01])
        elif func == 'vrg_Z':
            fitfunc = f_vrg
            bounds  = ([-.1,line-tol_aa,0. ,0,  1,-1.3,0,-.01], # y=-.01 = larger EWs
                       [ .0,line+tol_aa,1.5,1,410, 0. ,2, .01])
        



        #============================= Line fitting ============================
        i = 0; width_i = width
        p0_considered_def = []
        while i < iter:
            # Extracting the window of the spectrum
            if type(line_c) == float: window = (self.wave >= line-width_i/2) & (self.wave <= line+width_i/2)
            else:  window = (self.wave >= min(line_c)-width_i/2) & (self.wave <= max(line_c)+width_i/2)
            #===================== Continuum determination ====================
                
            if not any(window):
                print('Line %sA not in spectra.\n' % line)
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
                    ax.plot(wave,np.ones(len(wave)))
                    c_fit = np.poly1d(np.polyfit(wave, flux, 1))
                    continuum_i = c_fit(wave)
                    # Final normalization of the iteration
                    if normalization is True:
                        flux_norm = flux / continuum_i
                    else: 
                        flux_norm = flux
                    ax.plot(wave, flux_norm, c='b', lw=.5)
                    ax.set_title(self.id_star, fontsize=8)

                    #ax.set_yticks([])
                    ax.set_xlabel('$\lambda$ $[\AA]$', size=13)
                    ax.set_ylabel('Normalized flux', size=13)
                    ax.tick_params(direction='in', top='on', labelright=True, labelleft=False)
                    ax.figure.subplots_adjust(top=.9, bottom=.12, right=.88, left=.08)
                plt.show(block=False)                
                return fitsol

            flux = self.flux[window]
            wave = self.wave[window]
            if fit_range: 
                flux = flux[(wave >= fit_range[0]) & (wave <= fit_range[1])]
                wave = wave[(wave >= fit_range[0]) & (wave <= fit_range[1])] 
            dlam_mean = (wave[-1]-wave[0])/(len(wave)-1)
            # Auto-resampling
            #if dlam_mean >= 0.025 and not 'log' in star.filename:
            #    factor = dlam_mean/0.025; star.resamp(factor)

            # Find regions of the continuum to use during normalization (4 iter)
            '''
            for j in range(4):
                if j == 0:
                    mask_i = ~np.isnan(flux)
                else:
                    mask_i = ~sigma_clip(flux/continuum_i, maxiters=None,
                        sigma_lower=1.4, sigma_upper=2.5, axis=-1).mask #1.4 before

                c_fit = np.poly1d(np.polyfit(wave[mask_i], flux[mask_i], 1))
                continuum_i = c_fit(wave)
            '''
            # Final normalization of the iteration
            if normalization is True:         
                c_fit = np.poly1d(np.polyfit(wave, flux, 1))
                continuum_i = c_fit(wave)
                flux/continuum_i
                mask_i = ~sigma_clip(flux, maxiters=None,
                   sigma_lower=sigma_tres, sigma_upper=sigma_tres, axis=-1).mask #1.4 before  
                
                c_fit = np.poly1d(np.polyfit(wave[mask_i], flux[mask_i], 1))
                continuum_i = c_fit(wave)
                flux = flux/continuum_i
                flx,wv = flux[mask_i],wave[mask_i]
    
                while len(flx[~sigma_clip(flx, maxiters=None,sigma_lower=sigma_tres, sigma_upper=sigma_tres, axis=-1).mask])<len(flx):
                    
                    mask_i = ~sigma_clip(flx, maxiters=None,sigma_lower=sigma_tres, sigma_upper=sigma_tres, axis=-1).mask
                    flx,wv = flx[mask_i],wv[mask_i]
                    c_fit = np.poly1d(np.polyfit(wv, flx, 1))
                    continuum_i = c_fit(wave)
                    flux = flux/continuum_i
                mask_i = np.zeros(len(wave), dtype=bool)
                mask_i[np.intersect1d(wave, wv, return_indices=True)[1]] = True
            else: mask_i = np.ones(len(wave), dtype=bool); continuum_i = np.ones(len(wave))
            mask = mask_i
            flux_norm_i = flux
            #========================= Fitting the line ========================ç
            
            try:    
                
                if func == 'g':
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
                    #elif func == 'gc' and p0 != None and type(p0[1]) == list: 
                    
                else: 
                    p0_considered = np.array(p0)
                    bounds_considered = np.array(bounds[0]), np.array(bounds[1])
                popt_i = curve_fit(fitfunc, wave, flux_norm_i, bounds=bounds_considered, p0=p0_considered)[0]
                print('result',popt_i)
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
                FWHM_min = np.max([3*dlam_mean,3/4*dlamb])
                print(FWHM_min)
                FWHM_min = min(FWHM_min,0.05)
                if FWHM_min < FWHM < FWHM_max:
                    flux_norm = flux_norm_i; continuum = continuum_i; mask = mask_i
                    flux_fit = flux_fit_i; popt = popt_i; width = width_i
                    i = i + 1; width_i = FWHM*7

                    # Possible centers 
                    closest_WaveCentre_A,closest_WaveCentre_lam0 = [],[]
                    EW_dict = {}
                    # We compute the gaussian corresponding to the studied line. 
                    if func == 'g':
                        flux_fit_EW = f_gaussian1(wave, *popt_i[:3])
                    elif func == 'gc':
                        flux_fit_EW = gaussian_convolved(wave,popt_i[0],popt_i[len(popt_i)//2-2],*popt_i[-4:])
                        if len(popt_i[len(popt_i)//2-2:-4])>1:
                            print(popt_i[len(popt_i)//2-2:-4])
                            closest_WaveCentre_A.append(popt_i[0])
                            closest_WaveCentre_lam0.append(popt_i[len(popt_i)//2-2])
                            for index in range(1,len(popt_i[len(popt_i)//2-2:-4])):
                                # if abs(popt_i[len(popt_i)//2-2+index]-popt_i[len(popt_i)//2-2])<0.4:
                                if abs(float(p0_considered[len(p0_considered)//2-2])\
                                       -float(p0_considered[len(p0_considered)//2-2+index]))<0.4:
                                    # If difference in the THEORETICAL central value is lower than 0.4; basically the same line
                                    closest_WaveCentre_A.append(popt_i[0+index])
                                    closest_WaveCentre_lam0.append(popt_i[len(popt_i)//2-2+index])
                            flux_fit_EW = gaussian_convolved_general(wave,\
                            *closest_WaveCentre_A,*closest_WaveCentre_lam0,*popt_i[-4:])
                            print('Double line considered: ',closest_WaveCentre_A,closest_WaveCentre_lam0)
                            for critical_line in '4510.965','4514.854','4518.143','4523.561':
                                if (critical_line in p0_considered and \
                                    abs(float(p0_considered[len(p0_considered)//2-2])-float(critical_line))>1e-4) or\
                                    (round(float(critical_line),3) in p0_considered and \
                                    abs(float(p0_considered[len(p0_considered)//2-2])-float(critical_line))>1e-4):
                                    #popt_i = popt_i.astype(float)
                                    print('We are in',critical_line)
                                    cw,idx = fm.find_nearest(popt_i, float(critical_line))
                                    flux_fit_EW_C = gaussian_convolved(wave,popt_i[idx-(len(popt_i)//2-2)],cw,*popt_i[-4:])
                                    print('CHECK 1')
                                    if len(popt_i[len(popt_i)//2-2:-4])>1:
                                        closest_WaveCentre_A_C,closest_WaveCentre_lam0_C = [],[]
                                        for index in range(1,len(popt_i[len(popt_i)//2-2:-4])):
                                            # if abs(popt_i[len(popt_i)//2-2+index]-popt_i[len(popt_i)//2-2])<0.4:
                                            if abs(float(critical_line)-float(p0_considered[len(p0_considered)//2-2+index]))<0.4:
                                                print('In if')
                                                # If difference in the THEORETICAL central value is lower than 0.4
                                                closest_WaveCentre_A_C.append(popt_i[0+index])
                                                closest_WaveCentre_lam0_C.append(popt_i[len(popt_i)//2-2+index])
                                        print('CHECK 2')
                                        flux_fit_EW_C = gaussian_convolved_general(wave,\
                                        *closest_WaveCentre_A_C,*closest_WaveCentre_lam0_C,*popt_i[-4:])
                                        flux_fit_EW_C_min = gaussian_convolved_general(wave,\
                                        *closest_WaveCentre_A_C*(1+3/self.snr),*closest_WaveCentre_lam0_C,*popt_i[-4:])
                                        flux_fit_EW_C_max = gaussian_convolved_general(wave,\
                                        *closest_WaveCentre_A_C*(1+3/self.snr),*closest_WaveCentre_lam0_C,*popt_i[-4:])
                                        print(*closest_WaveCentre_A_C,*closest_WaveCentre_lam0_C)
                                    ####### Computing critical line EW  #######
                                    EW = 0.5*abs(fsum((wave[wl-1] - wave[wl])*((1 - flux_fit_EW_C[wl-1]) +
                                         (1 - flux_fit_EW_C[wl])) for wl in range(1, len(flux_fit_EW_C))))
                                    EW = round(1000*EW)
                                    EW_dict[float(critical_line)] = EW
                            '''
                            Three lines which only fits sometimes: 4510.965,4514.854,4523.561
                            '''
                elif FWHM < FWHM_min:
                    print('WARNING: FWHM(%.1f) < minimum FWHM for %.3fA' % (FWHM,line))
                    break
                elif FWHM > FWHM_max:
                    print('WARNING: FWHM(%.1f) > maximum FWHM for %.3fA ' % (FWHM,line))
                    break
                
            except: 
                print('Not feasable')
                break
        fitsol['lines'] = p0_considered_def
        print('Considered lines for ',line,fitsol['lines'])
        
        #======================== Checking final results =======================
        #window = (self.wave >= line-width/2.) & (self.wave <= line+width/2.)
        
        if type(line_c) == float: window2 = (self.wave >= line-width/2) & (self.wave <= line+width/2)
        else:  window2 = (self.wave >= min(line_c)-width/2) & (self.wave <= max(line_c)+width/2)
    
        flux = self.flux[window2]
        wave = self.wave[window2]
        if fit_range: 
            flux = flux[(wave >= fit_range[0]) & (wave <= fit_range[1])]
            wave = wave[(wave >= fit_range[0]) & (wave <= fit_range[1])] 
        
        if len(self.wave[window]) != len(wave) and normalization is True:  
            print('True. Recomputing...')
            c_fit = np.poly1d(np.polyfit(wave, flux, 1))
            continuum_i = c_fit(wave)
            flux/continuum_i
            mask_i = ~sigma_clip(flux, maxiters=None,
               sigma_lower=sigma_tres, sigma_upper=sigma_tres, axis=-1).mask #1.4 before  
            
            c_fit = np.poly1d(np.polyfit(wave[mask_i], flux[mask_i], 1))
            continuum_i = c_fit(wave)
            flux = flux/continuum_i
            flx,wv = flux[mask_i],wave[mask_i]
    
            while len(flx[~sigma_clip(flx, maxiters=None,sigma_lower=sigma_tres, sigma_upper=sigma_tres, axis=-1).mask])<len(flx):
                
                mask_i = ~sigma_clip(flx, maxiters=None,sigma_lower=sigma_tres, sigma_upper=sigma_tres, axis=-1).mask
                flx,wv = flx[mask_i],wv[mask_i]
                c_fit = np.poly1d(np.polyfit(wv, flx, 1))
                continuum_i = c_fit(wave)
                flux = flux/continuum_i
            mask_i = np.zeros(len(wave), dtype=bool)
            mask_i[np.intersect1d(wave, wv, return_indices=True)[1]] = True
        else: mask_i = np.ones(len(wave), dtype=bool); continuum_i = np.ones(len(wave))
        mask = mask_i

        if i == 0:
            if info is True:
                print('Problem in spectrum %s' % self.filename)
                print('Line %sA could not be fitted or does not exist.\n' % line)
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

                ax.plot(wave,np.ones(len(wave)))
                c_fit = np.poly1d(np.polyfit(wave, flux, 1))
                continuum_i = c_fit(wave)
                # Final normalization of the iteration
                #if normalization is True:
                #    flux_norm = flux / continuum_i
                #else: flux_norm = flux 
                flux_norm = flux 
                ax.plot(wave, flux_norm, c='b', lw=.5)
                ax.set_title(self.id_star, fontsize=8)
                ax.plot(wave, np.where(mask==False, 1, np.nan) + 0.005, 'k', lw=.5)


                #ax.set_yticks([])
                ax.set_xlabel('$\lambda$ $[\AA]$', size=13)
                ax.set_ylabel('Normalized flux', size=13)
                ax.tick_params(direction='in', top='on', labelright=True, labelleft=False)
                ax.figure.subplots_adjust(top=.9, bottom=.12, right=.88, left=.08)
            plt.show(block=False)
            return fitsol

        if popt_i[0]<0:
            #line_f = wave[flux_fit == min(flux_fit)][0]
            line_f = fm.find_nearest(popt_i, float(line))[0]
        else:
            #line_f = wave[flux_fit == max(flux_fit)][0]
            line_f = fm.find_nearest(popt_i, float(line))[0]
        
        if abs(line - line_f) > tol_aa:
            if info is True:
                print('Line %sA found outside tolerance.\n' % line)
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

                ax.plot(wave,np.ones(len(wave)))
                c_fit = np.poly1d(np.polyfit(wave, flux, 1))
                continuum_i = c_fit(wave)
                # Final normalization of the iteration
                if normalization is True:
                    flux_norm = flux / continuum_i
                else: flux_norm = flux
                ax.plot(wave, flux_norm, c='b', lw=.5)
                ax.set_title(self.id_star, fontsize=8)
                ax.plot(wave, np.where(mask==False, 1, np.nan) + 0.005, 'k', lw=.5)

                #ax.set_yticks([])
                ax.set_xlabel('$\lambda$ $[\AA]$', size=13)
                ax.set_ylabel('Normalized flux', size=13)
                ax.tick_params(direction='in', top='on', labelright=True, labelleft=False)
                ax.figure.subplots_adjust(top=.9, bottom=.12, right=.88, left=.08)
            plt.show(block=False)
            
            return fitsol

        RV_A   = round((line_f - line), 3)
        RV_kms = round(((line_f - line)/line)*cte.c/1000, 2) # max precision is 100 m/s
        line_f = round(line_f, 3)
        if info is True:
            print('Line %sA found at %.3fA -> RV: %.1fkm/s\n' % (line,line_f,RV_kms))

        #=========================== Calculate the EW ==========================
        # stackoverflow.com/questions/34075111/calculate-equivalent-width-using-python-code
        # When emission is considered abs should be removed and 1-flux_fit -> flux_fit-1
        # Compute the EW from the flux of the gaussian FROM the studied line. 
        EW = 0.5*abs(fsum((wave[wl-1] - wave[wl])*((1 - flux_fit_EW[wl-1]) +
             (1 - flux_fit_EW[wl])) for wl in range(1, len(flux_fit_EW))))  
        EW = round(1000*EW)
        EW_dict[line_f] = EW
        print('EW='+str(EW))
        print('EW for ',line_f,EW_dict)
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

        #===================== Calculate the SNR continuum =====================
        print(len(flux_norm),len(mask),wave[0],wave[-1])
        sigma_cont = np.std(flux_norm[mask])
        snr = int(1/sigma_cont)

        #======= Calculate the SNR continuum with only removing the line =======
        mask2 = (wave < line_f - FWHM) + (wave > line_f + FWHM)
        sigma_cont2 = np.std(flux_norm[mask2])
        
        if len(flux_norm[mask2])<=1: snr2 = np.nan
        else: snr2 = int(1/sigma_cont2)
        
        #============================= Quality value ===========================
        q_fit = 1/np.std(flux_norm[flux_fit<(1-.2*depth)]/flux_fit[flux_fit<(1-.2*depth)]) #simple
        q_fit = round(q_fit, 3)
        
        #========================== EW error stimation ========================       
        if EW_err_k:
            continuum_ex,EW_MCMC = [],[]
            start = time.time()
            while i<iter_err:
                try:
                    cont_err = np.random.normal(loc = 1,scale=sigma_cont)
                    popt_err = curve_fit(fitfunc, wave, flux_norm_i/cont_err, \
                                    bounds=bounds_considered, p0=p0_considered)[0]
                    if func == 'g': flux_fit_EW = f_gaussian1(wave, *popt_err[:3])
                    elif func == 'gc': 
                        flux_fit_EW = gaussian_convolved_general(wave,\
                                    *closest_WaveCentre_A,*closest_WaveCentre_lam0,*popt_err[-4:])
                    EW_err = 0.5*abs(fsum((wave[wl-1] - wave[wl])*((1 - flux_fit_EW[wl-1]) +
                         (1 - flux_fit_EW[wl])) for wl in range(1, len(flux_fit_EW))))
                    EW_err = round(1000*EW_err)      
                    EW_MCMC.append(EW_err)
                    i += 1     
                    #print(i)
                except: pass
            print('Continuum error is ',sigma_cont)
            print('EW error=',np.std(EW_MCMC))
            fig,ax = plt.subplots()
            ax.hist(EW_MCMC,bins=10,density=True)
            ax.grid()
            print('final time:',round(time.time()-start,2))
            ### Limit values error stimation : ###
            popt_err = curve_fit(fitfunc, wave, flux_norm_i/(1+sigma_cont), bounds=bounds_considered, p0=p0_considered)[0]
            if func == 'g': flux_fit_EW = f_gaussian1(wave, *popt_err[:3])
            elif func == 'gc': 
                flux_fit_EW = gaussian_convolved_general(wave,\
                            *closest_WaveCentre_A,*closest_WaveCentre_lam0,*popt_err[-4:])
            EW_inf = 0.5*abs(fsum((wave[wl-1] - wave[wl])*((1 - flux_fit_EW[wl-1]) +
                 (1 - flux_fit_EW[wl])) for wl in range(1, len(flux_fit_EW))))
            EW_inf = round(1000*EW_inf)             
            popt_err = curve_fit(fitfunc, wave, flux_norm_i/(1-sigma_cont), bounds=bounds_considered, p0=p0_considered)[0]
            if func == 'g': flux_fit_EW = f_gaussian1(wave, *popt_err[:3])
            elif func == 'gc': 
                flux_fit_EW = gaussian_convolved_general(wave,\
                            *closest_WaveCentre_A,*closest_WaveCentre_lam0,*popt_err[-4:])
            
            EW_max = 0.5*abs(fsum((wave[wl-1] - wave[wl])*((1 - flux_fit_EW[wl-1]) +
                 (1 - flux_fit_EW[wl])) for wl in range(1, len(flux_fit_EW))))
            EW_max = round(1000*EW_max)          
            print('EW limit error=',(EW_inf-EW_max)/2)
            ax.axvline(EW_inf,c='r'),ax.axvline(EW_max,c='r')
            ax.set_xlabel(r'EW [m$\AA$]')
        #================================ Plot =================================
        rest = flux_fit - flux_norm
        dEW = EW * 3*(np.std(rest) + abs(np.mean(rest)))
        #print(dEW)
        #plt.figure()
        #plt.plot(rest)
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
            ax.plot(wave, continuum, c='r', lw=.5)
            ax.plot(wave, flux_norm, c='b', lw=.5)
            ax.plot(wave, np.where(mask==False, 1, np.nan) + 0.005, 'k', lw=.5)
            #ax.plot(wave, np.where(mask2==False, 1, np.nan) + 0.005, 'grey', lw=.5)
            if func=='gc':
                ax.set_title('%s | %.2f | RV: %d | EW: %d | dEW: %.2f | FWHM: %.2f | SNR: %d | vsini: %.2f | vmac: %.2f' %
                   (self.id_star,line_f,RV_kms,EW,dEW,FWHM,snr,popt_i[-3],popt_i[-1]), fontsize=8)               
            else:
                ax.set_title('%s | %.2f | RV: %d | EW: %d | dEW: %.2f | FWHM: %.2f | FW34-14: %.2f | SNR: %d' %
                    (self.id_star,line_f,RV_kms,EW,dEW,FWHM,FW34_14,snr), fontsize=8)
            
            
            #ax.set_yticks([])
            ax.set_xlabel('$\lambda$ $[\AA]$', size=13)
            ax.set_ylabel('Normalized flux', size=13)
            ax.tick_params(direction='in', top='on', labelright=True, labelleft=False)
            ax.figure.subplots_adjust(top=.9, bottom=.12, right=.88, left=.08)
            ax.set_xlim(max(line-3*FWHM,plt.xlim()[0]),min(line+3*FWHM,plt.xlim()[1]))
            
        plt.show(block=False)

        #=======================================================================
        #================= Packing the results in a dictionary =================
        fitsol = {'sol':1, 'line':line_f, 'RV_A':RV_A, 'RV_kms':RV_kms,
                   'EW':EW, 'EW_err':dEW, 'FWHM':FWHM, 'FW34_14':FW34_14, 'depth':depth,
                   'q_fit':q_fit, 'snr':snr, 'snr2':snr2, 'lines': popt_i,#p0_considered_def,
                   'EW_dict':EW_dict, 'line_EW_studied':line}
        if EW_err_k:
            fitsol['EW_err'] = round(np.std(EW_MCMC),2)
        if info is True:
            print(line,popt_i)
        if func == 'gc' and type(line_c) != float:
            popt = popt[:len(line_c)],popt[len(line_c):2*len(line_c)],*popt[-4:]
        for f_par,par in zip(fit_dic[func.split('_')[0]], popt):
            fitsol[f_par] = np.around(par, 3)

        if outfit == True:
            fitsol['wave'] = wave
            #fitsol['flux'] = flux
            fitsol['flux_norm'] = flux_norm
            fitsol['flux_fit'] = flux_fit
        return fitsol

        # Theorerical FWHM:
        #if   func == 'g': jFWHM = 2*np.sqrt(2*np.log(2))*popt[2]
        #elif func == 'l': jFWHM = 2*abs(popt[2])
        #elif func == 'v': jFWHM = 2*(.5346*popt[3]+np.sqrt(.2166*(popt[3]**2)+popt[2]**2))
        #elif func == 'r': jFWHM = 1.7*popt[3]*line*1000/cte.c
        #jFWHM = round(jFWHM, 3)


    def snrcalc(self, zone='B'):

        '''
        Function to calculate the Signal to Noise Ratio in different regions of the
        spectra if available.

        Parameters
        ----------
        zone : str, optional
            Select the zone to calculate the spectra.
                'b'/'B'     -> 4000-5000 A
                'v'/'V'     -> 5000-6000 A
                'r'/'R'     -> 6000-7000 A
                'all'/'ALL' -> 4000-7000 A

        Returns
        -------
        Measured signal-to-noise ratio value.
        '''

        if zone in ['b','B']:
            mask = (self.wave > 4000) & (self.wave < 5000)
        elif zone in ['v','V']:
            mask = (self.wave > 5000) & (self.wave < 6000)
        elif zone in ['r','R']:
            mask = (self.wave > 6000) & (self.wave < 7000)
        elif zone in ['all','ALL']:
            mask = (self.wave > 4000) & (self.wave < 7000)

        lambda0 = np.mean(self.wave[mask])
        resol = 10000

        sigma = lambda0/(2.35482*float(resol))

        gauss = f_gaussian(np.arange(-5*sigma, 5*sigma, self.dlam), sigma)
        kernel = gauss/np.trapz(gauss)

        convoluted = 1 + convolve(self.flux[mask] - 1, kernel, mode='same')

        flux_norm = self.flux[mask]/convoluted

        snr_all = []
        for gap in findlist('snr_gaps.txt'):
            lwl,rwl = [float(i) for i in gap.split('-')]

            flux_norm_i = flux_norm[(self.wave[mask] >= lwl) & (self.wave[mask] <= rwl)]

            if len(flux_norm_i) == 0:
                continue

            sig_clip = 3
            std = np.std(flux_norm_i)
            flux_clean = np.where(abs(flux_norm_i - 1) > sig_clip*std, np.nan, flux_norm_i)

            snr_all.append(1/np.nanstd(flux_clean))

        self.snr = np.nanmean(snr_all)

        if not np.isnan(np.nanmean(snr_all)):
            return int(round(np.nanmean(snr_all)))


    def cosmic(self, method='zscore',
        dmin=0.05, zs_cut=5,
        wl_split=100, ker_sig=2, ker_iter=3, sig_g=None,
        protect_em_lines=True):

        '''
        Function to remove cosmic rays in the spectra by different approaches.

        Parameters
        ----------
        method : str, optional
            Method for the cosmic ray removal strategy. Only zscore (def) or kernel.

        dmin : int/float, optional
            Minium distance between flux and cleaned flux to consider for replacement.
            Default is 0.05.

        zs_cut : int/float, optional
            Threshold value used in the zscore method for finding rays. Default is 4.
            Tip: For noisy spectra, rise this value up to 7-9.

        wl_split : int/float, optional
            In the 'kernel' method, wavelenght size used to split the spectrum before
            applying the cosmic removal. Default is 100.

        ker_sig : float, optional
            Sigma clipping value used to remove rays. Default is 2.

        ker_iter : int, optional
            Number of iterations of the sigma clipping to remove cosmic rays in the kernel method.
            Default is 3.

        sig_g : float, optional
            Sigma of the gaussian function used to construct the kernel.
            Default is the theoretical sigma based on wavelenght and resolution.

        protect_em_lines : boolean, optional
            If True, some emission lines will be masked from cosmic rays removal.
            Default is True.

        Returns
        -------
        Nothing, but the flux is replaced and cleaned from rays.
        '''

        if method == 'zscore':
            # www.towardsdatascience.com/removing-spikes-from-raman-spectra-8a9fdda0ac22

            # First we calculated (nabla)x(i):
            delta_flux = np.diff(self.flux)
            median_int = np.median(delta_flux)
            mad_int = np.median([np.abs(delta_flux - median_int)])
            modified_z_scores = 0.6745 * (delta_flux - median_int) / mad_int
            # The multiplier 0.6745 is the 0.75th quartile of the standard normal
            # distribution, to which the median absolute deviation converges to.
            modified_z_scores =  np.concatenate(([0], np.abs(modified_z_scores)))

            flux_clean = np.where(modified_z_scores > zs_cut, np.nan, self.flux)

        elif method == 'kernel':

            wl_split = 100 # Range in angstroms in which the spectrum will be initially splitted

            wl_range = self.wave[-1]-self.wave[0]
            if wl_range < wl_split:
                wl_split = wl_range

            if 1 < wl_range/wl_split <= 2:
                wl_split = wl_range/2 + 1
            elif wl_range/wl_split > 2:
                wl_split += (wl_range % wl_split) / int(wl_range/wl_split)

            n = int(wl_range/wl_split)

            flux_clean = []
            for i in range(n):

                mask = (self.wave >= self.wave[0]+i*wl_split) & (self.wave < self.wave[0]+(i+1)*wl_split)
                if i == range(n)[-1]:
                    mask[-1] = True # To catch the last flux value from not using <= above

                resolution = float(self.resolution) # resolution = 5000 # Empirical optimal value
                dlam = self.dlam # dlam = 0.2564975 # Empirical optimal value

                if sig_g is None:
                    lambda0 = np.mean(self.wave[mask])
                    sig_g = lambda0/(2.35482*resolution)
                else:
                    sig_g = float(sig_g)

                x = np.arange(-5*sig_g, 5*sig_g+dlam, dlam)
                gauss = f_gaussian(x,sig_g)
                kernel = gauss/np.trapz(gauss)

                convoluted = 1 + convolve(self.flux[mask] - 1, kernel, mode='same')

                flux_norm = self.flux[mask]/convoluted

                for i in range(ker_iter):
                    std = np.nanstd(flux_norm)
                    flux_norm = np.where(abs(flux_norm - 1) > ker_sig*std, np.nan, flux_norm)

                flux_clean = np.concatenate([flux_clean,np.where(np.isnan(flux_norm), np.nan, self.flux[mask])])

        nans = np.isnan(flux_clean)
        x = lambda z: z.nonzero()[0]
        flux_clean[nans] = np.interp(x(nans), x(~nans), flux_clean[~nans])

        # Recover the original flux in the regions of the spectrum where telluric lines are
        flux_clean = np.where(
              ((self.wave > 3932.5) & (self.wave < 3934.5)) # Not sure if this is exactly telluric line
            | ((self.wave > 5885.0) & (self.wave < 5900.0))
            | ((self.wave > 6865.0) & (self.wave < 7035.0))
            | ((self.wave > 7160.0) & (self.wave < 7340.0))
            | ((self.wave > 7585.0) & (self.wave < 7700.0))
            | ((self.wave > 8116.0) & (self.wave < 8380.0))
            | ((self.wave > 8915.0) & (self.wave < 9220.0))
               ,self.flux, flux_clean)

        # Recover the original flux in the regions of the spectrum where emission lines are
        if protect_em_lines == True:
            for wl_em in [3967.79,4958.911,5006.843,6300.304,6548.04,6583.46,6716.44,6730.82]:
                mask = (self.wave > wl_em-0.8) & (self.wave < wl_em+0.8)
                flux_clean[mask] = self.flux[mask]

        # Recover the original flux if the difference is lower than dmin value
        flux_clean = np.where(
            (self.flux > flux_clean) & (self.flux - flux_clean > dmin), flux_clean, self.flux)

        self.flux = flux_clean

        return None


    def cosmetic():

        '''
        IN DEVELOPMENT - Function to remove cosmetic defects recurrently found in the
        spectra.
        '''

        # Correct from FEROS issues for window in
        #[[4506,4507.5],[4693.7,4696],[4795.,4797.],[4900.7,4902.3],[6636.3,6638.3]]
        return None


    def degrade(self, resol, profile='g', vsini=None, vmac=None):

        '''
        Function to degrade a spectrum to a certain resolution by convolving it to a
        gaussian (pure degradation) or to account for rotational+macroturbulence effect
        for example if a synthetic spectrum is loaded.

        Parameters
        ----------
        resol : int/float, optional
            Resolution of the gaussian profile used to degrade the spectrum.

        profile : str
            Use 'g' for gaussian profile convolution (Default).
            Use 'rotmac' for rotational+macroturbulence profile convolution.

        vsini : int/float, optiomal
            Value of vsini. Only valid for rotational+macroturbulence profile.

        vmac : int/float, optiomal
            Value of vmac. Only valid for rotational+macroturbulence profile.

        Returns
        -------
        Nothing, but the flux is replaced by the degraded one.
        '''

        lambda0 = np.mean(self.wave)

        if profile == 'g' and (vsini==None and vmac==None):
            sigma = lambda0/(2.35482*float(resol))

            x = np.arange(-10*sigma, 10*sigma+self.dlam, self.dlam)
            gauss = f_gaussian(x,sigma)
            kernel = gauss/np.trapz(gauss)
            self.resolution = resol

        elif profile == 'rotmac' and (vsini!=None and vmac!=None):
            x = np.arange(-9, 9+self.dlam, self.dlam)
            rotmac = f_rotmac(x, lambda0, vsini, vmac)
            kernel = rotmac/np.trapz(rotmac)

        convoluted = 1 + convolve(self.flux - 1, kernel, mode='same')

        self.flux = convoluted


    def resamp(self, dlam, lwl=None, rwl=None, method='linear'):

        '''
        Function to resample a spectrum into a fixed delta-lambda and wavelenght range.

        Parameters
        ----------
        dlam : float/int
            New delta lambda to be used for the output spectra.

        lwl : float/int, optional
            Enter the forced initial wavelenght to be used during interpolation.
            If None, the original initial wavelenght will be used.

        rwl : float/int, optional
            Enter the forced final wavelenght to be used during interpolation.
            If None, the original final wavelenght will be used.

        method : str, optional
            Enter the interpolation method to be used. See doc for np.interp1d.
            Default is 'linear'.

        Returns
        -------
        Nothing, but the spectrum (wavelenght,flux) is resampled.
        '''

        try:
            float(dlam)
        except:
            print('Input should be float or integrer.'); return None

        self.dlam = dlam

        if dlam > np.mean(self.wave)/self.resolution/3:
            # It is divided by 3 to at least have 3 pixels in a gaussian
            print('WARNING: The new delta lambda implies lossing information...')

        if lwl is None or lwl < self.wave[0]:
            lwl = self.wave[0]
        if rwl is None or rwl > self.wave[-1]:
            rwl = self.wave[-1]

        f = interp1d(self.wave, self.flux, kind=method, fill_value='extrapolate')
        self.wave = np.arange(lwl, rwl+self.dlam, self.dlam)
        self.flux = f(self.wave)


    def export(self, tail='', extension='.ascii'):

        '''
        Function to export the current wavelenght and flux of the spectrum in the class
        into an ascii file.

        Parameters
        ----------
        tail : str, optional
            Tail of the file added before the extension for its identification.
            Default is ''.

        extension : str, optional
            Extenstion of the output file. Default is '.ascii'.

        Returns
        -------
        Nothing, but the ascii file is exported.
        '''

        filename = self.filename.replace('.fits', '')
        np.savetxt(maindir+'tmp/%s' % (filename + tail + extension),
            np.c_[self.wave,self.flux], fmt=('%.4f','%.6f'),
            header='lambda    flux', comments='')

        return None


    def plotline(self, lines, width=10, ylim=None):

        '''
        Function to create a plot around a spectral line or lines.

        Parameters
        ----------
        lines : float, str
            Enter the wavelenght(s) of the line(s) to plot, either in a coma-separated
            string, or in a .txt/.lst file containing the lines.

        width : int, optional
            Sets the width in [AA] where the line fits well in. Default is 10.

        ylim : tuple/list, optional
            Sets the y-limits for the plot. Input must be like "[ymin,ymax]".

        Returns
        -------
        Nothing, but the plots are generated.
        '''

        self.get_spc()

        lines,elements,_ = findlines(lines)
        if len(lines) > 1:
            nrows = ncols = int(round(np.sqrt(len(lines)), 0))
        else:
            nrows = ncols = 1

        for line,element,nplot in zip(lines, elements, range(len(lines))):

            mask = (self.wave > line - width/2) & (self.wave < line + width/2)

            if len(lines) > 1:
                plt.subplot(nrows, ncols, nplot + 1)
                plt.xticks([round(line - width/3, 1),round(line, 1),round(line + width/3, 1)])
                plt.title(element, fontsize=6, pad=1)

            plt.plot(self.wave[mask], self.flux[mask], lw=.3, label=self.id_star+' '+self.SpC)
            plt.tick_params(direction='in', top='on')

            if ylim is not None and (type(ylim) is list or type(ylim) is tuple):
                plt.ylim(ylim)

            if len(lines) == 1:
                plt.xlabel('$\lambda$ $[\AA]$', size=13)
                plt.ylabel('Normalized flux', size=13)

            plt.tight_layout()

        plt.legend()
        plt.show(block=False)

        return None


    def plotspec(self, lwl=3800, rwl=8000, lines=None, ylim=None, path='sp_lines/'):

        '''
        Function to create a plot of a portion of the spectra and optionally overplot
        tabulated spectral lines in that range taken from a database.

        Parameters
        ----------
        lwl : float, optional
            Sets the start wavelenght of the spectrum.

        rwl : float, optional
            Sets the end wavelenght of the spectrum.

        lines : str, optional
            Use to overplot position of spectral lines. Current options are:
            ALL, ALLOB, 35-10K, 35K, 30K, 25K, 20K, 15K, 10K

        ylim : tuple/list, optional
            Sets the y-limits for the plot. Input must be like "[ymin,ymax]".

        path : str, optional
            Path to where the spectral line files are located. Default inside sp_lines/.

        Returns
        -------
        Nothing, but the plots are generated.
        '''

        self.get_spc()

        if lwl < min(self.wave):
            lwl = min(self.wave)
        if rwl > max(self.wave):
            rwl = max(self.wave)

        mask = (self.wave > lwl) & (self.wave < rwl)

        if lines != None:

            try:
                depth = max(self.flux[mask]) - min(self.flux[mask]) # or 1-min
            except:
                print('Problem finding max/min in masked flux.')
                return None

            if self.rv0 == 0:
                print('Spectrum not corrected from RV, lines will have offset.')

            if  lines == 'ALL':
                table = findtable('ALL_all.txt', delimiter=',', path=path)
                table = table[table['-lg(gf)'] > -1]
            elif lines == 'ALLOB':
                table = findtable('ALL_OBs_n4+.txt', delimiter=',', path=path)
                table = table[table['-lg(gf)'] > -1]
            elif lines in ['35-10K','35K','30K','25K','20K','15K','10K']:
                table = findtable('%s.fits' % lines, path=path)
                # www.lsw.uni-heidelberg.de/projects/hot-stars/websynspec.php

            table = table[(table['wl_air'] >= lwl) & (table['wl_air'] <= rwl)]

            if '-lg(gf)' in table.columns:
                table['width'] = 10**table['-lg(gf)']/np.max(10**table['-lg(gf)']) * 3
                # 10**gf/5 empiric way to draw thicker lines for instense lines
            elif 'strength' in table.columns:
                table['strength'] = [i if i<800 else 800 for i in table['strength']]
                table['width'] = table['strength']/np.max(table['strength']) * 3

            at_color = {'HI':'gray', 'HeI':'turquoise', 'OI':'r', 'NI':'b', 'CI':'k', 'SI':'gold',
                'Si':'tan', 'Mg':'g', 'Fe':'chocolate', 'Ne':'teal', 'Al':'rosybrown'}

            for line in table:
                try: c = at_color[line['spc'].replace(' ','')[:2]]
                except: c = 'dimgray'

                plt.plot([line['wl_air'],line['wl_air']], [1.008-depth,np.median(self.flux[mask])],
                c=c, linestyle='dotted', lw=line['width'])

                # depth line mask = depth deepest line
                plt.text(line['wl_air'],1.004-depth, line['spc'], c=c, size=6, rotation=-90, clip_on=True)

        plt.plot(self.wave[mask], self.flux[mask], lw=.3, label=self.id_star+' '+self.SpC)
        plt.tick_params(direction='in', top='on')

        if ylim is not None and (type(ylim) is list or type(ylim) is tuple):
            plt.ylim(ylim)

        plt.xlabel('$\lambda$ $[\AA]$', size=13)
        plt.ylabel('Normalized flux', size=13)

        plt.legend()
        plt.tight_layout()
        plt.show(block=False)

        return None


# It now follows the functions describing the different fitting profiles:

def f_gaussian(x, sigma):
    return np.exp(-(x/sigma)**2/2)

def f_gaussian1(x, A, lam0, sigma):
    # A -> Amplitude;  lam0 -> center
    return A*np.exp(-(x - lam0)**2/(2*sigma**2)) + 1

def gaussian_convolved(x,A,lam0,sigma,vsini,beta,vmac,**kwargs):
    gaussian = f_gaussian1(x, A, lam0, sigma)
    broad_profile = Broadening_profile(x,lam0,vsini,beta,vmac)
    final_profile = broad_profile.convolution(gaussian-1)
    return final_profile[1] + 1

def gaussian_convolved_general(x,*params,**kwargs):
    f = 0 # All the three use the same sigma,vsini,beta,vmac
    for i in range(int(len(params)//2-2)):
        
        f += gaussian_convolved(x,params[i],params[int(len(params)//2-2)+i],params[-4],params[-3],params[-2],params[-1],**kwargs)-1
    #f += gaussian_convolved(x,A1,lams_01,sigma,vsini,beta,vmac,**kwargs)-1
    #f += gaussian_convolved(x,A2,lams_02,sigma,vsini,beta,vmac,**kwargs)-1
    f += 1
    return f    

def f_gaussian_general(x, *params):
    f = 0
    for i in range(len(params)//3):
        A,lam0,sigma = params[3*i],params[3*i+1],params[3*i+2]
        f += f_gaussian1(x, A, lam0, sigma)-1
    f += 1
    return f

def f_lorentzian(x, A, lam0, gamma, y):
    return A*gamma**2/((x - lam0)**2 + gamma**2) + y

def f_voigt(x, A, lam0, sigma, gamma, y):
    # sigma -> gaussian width; gamma -> lorentzian width
    # sigma = alpha / sqrt(2 * np.log(2))
    return A*np.real(wofz((x - lam0 + 1j*gamma)/sigma/np.sqrt(2)))/sigma/np.sqrt(2*np.pi) + y

def f_rot(x, A, lam0, sigma, vsini):
    G = A*np.exp(-(x - lam0)**2/(2*sigma**2))

    # Default value: beta=1.5 (epsilon=0.6) beta=epsilon/(1 - epsilon)
    eps = 0.6
    delta = 1000*lam0*vsini/cte.c
    doppl = 1 - ((x - lam0)/delta)**2

    R = A*(2*(1 - eps)*np.sqrt(doppl) + np.pi*eps/2.*doppl)/(np.pi*delta*(1 - eps/3))
    R = np.nan_to_num(R)

    return 1-convolve(G, R, mode='same')

def f_voigtrot(x, A, lam0, sigma, gamma, vsini, y):
    V = A*np.real(wofz((x-lam0+1j*gamma)/sigma/np.sqrt(2)))/sigma/np.sqrt(2*np.pi) + y

    eps = 0.6
    delta = 1000*lam0*vsini/cte.c
    doppl = 1 - ((x - lam0)/delta)**2

    R = A*(2*(1 - eps)*np.sqrt(doppl) + np.pi*eps/2.*doppl)/(np.pi*delta*(1 - eps/3))
    R = np.nan_to_num(R)

    return 1-convolve(V, R, mode='same')

def f_vrg(x, A, lam0, sigma, gamma, vsini, A2, sigma2, y):
    VG = A*np.real(wofz((x - lam0 + 1j*gamma)/sigma/np.sqrt(2)))/sigma/np.sqrt(2*np.pi) + y \
        + A2*np.exp(-(x - lam0)**2/(2*sigma2**2))

    eps = 0.6
    delta = 1000*lam0*vsini/cte.c
    doppl = 1 - ((x - lam0)/delta)**2

    R = A*(2*(1 - eps)*np.sqrt(doppl)+np.pi*eps/2.*doppl)/(np.pi*delta*(1 - eps/3))
    R = np.nan_to_num(R)

    return 1-convolve(VG, R, mode='same')

def f_rotmac(x, lam0, vsini=None, vmac=None):

    if vsini != None:
        # Rotational function:
        delta_R = 1000*lam0*vsini/cte.c
        doppl = 1 - (x/delta_R)**2

        eps = 0.6
        R = (2*(1 - eps)*np.sqrt(doppl) + np.pi*eps/2.*doppl)/(np.pi*delta_R*(1 - eps/3))
        R = np.nan_to_num(R)

        if vmac is None:
            return R

    if vmac != None:
        # Macroturbulence function:
        delta_M = 1000*lam0*vmac/cte.c
        A = 2/np.sqrt(np.pi)/delta_M

        x_2 = x[len(x)//2:]
        x_d = x_2/delta_M

        M_T = A*x_d*(-np.sqrt(np.pi)+np.exp(-x_d**2)/x_d+np.sqrt(np.pi)*erf(x_d))

        M = M_T # + M_R

        M = np.concatenate((M[::-1], M[1:]))

        if vsini is None:
            return M

    if vsini != None and vmac != None:
        return convolve(R, M, mode='same')
print('Read!')


#%% TEST
if __name__ == "__main__":
    theory,observations,fited_directly = False, False, False
    test = True
    
    
    from termcolor import colored    
    print(colored('RUNNING TEST', 'magenta', 'on_white',['bold']))
    import time
    start = time.time()
    plt.close('all')
    from termcolor import colored
    
    A = spec('HD151018')
    wave,flux,b = A.waveflux()
    
    
    line_c,vsin,vmac = 5410.8,67,65
    if theory:
        l0,vsini,vmac,R,beta,stp =5268.3,50,20,85000,1.5,0.005
        elems = 'C','N','O','H_He'
        elems_abundance = {'C': -3.7, 'N': -4.4, 'O': -3.25, 'H_He': -4.4}
        
        wave, flux = spectrum_generation('T340g410He10Q130b10CNO',elems,elems_abundance,'VT009')
        f1 = interp1d(wave,flux, kind='linear');wv = np.arange(start=wave[0],stop=wave[-1],step=0.05)
        flx = f1(wv) 
        print('Spectrum generation time:',time.time()-start)
        mask_orig = (wv>5220) & (wv<5310) 
        
        wvl1,flux1 = final_resolution(wv[mask_orig], flx[mask_orig]-1, l0, vsini, vmac, R, stp = np.mean(np.diff(wv)))
        mask_conv = (wvl1>5220) & (wvl1<5310)    
        
        print('New resolution',time.time()-start)
        #ej,class_ex = gaussian_convolved(wvl1[mask_conv],-0.5,l0,0.53,vsini,beta,vmac)
        
        popt_i = curve_fit(gaussian_convolved,wvl1[mask_conv],flux1[mask_conv]+1,\
        p0=[-0.17,l0,0.1,vsini,beta,vmac])[0]
        flux_fit_EW = gaussian_convolved(wvl1[mask_conv], *popt_i)#-0.17,l0,0.1,vsini,beta,vmac)#
        
        print(time.time()-start)
        
        EW_orig = sum((flx[mask_orig]-1)*np.mean(np.diff(wv[mask_orig])))*-1e3
        EW_convolved = sum(flux1[mask_conv]*np.mean(np.diff(wvl1[mask_conv])))*-1e3
        EW_fitted = sum(flux_fit_EW*np.mean(np.diff(wv[mask_orig])))*-1e3
        
        fitting_wave = np.linspace((wvl1[mask_conv][0]+wvl1[mask_conv][1])/2,\
                    (wvl1[mask_conv][-2]+wvl1[mask_conv][-1])/2,len(flux_fit_EW)-1)
        plt.figure()
        
        plt.plot(wv[mask_orig],flx[mask_orig],'k',label='FASTWIND model',alpha = 0.5)
        plt.plot(wvl1[mask_conv],flux1[mask_conv]+1,'b',label='FASTWIND convolved: vsini='+\
                  str(round(vsini,2))+'km/s, vmac='+str(round(vmac,2)))
        plt.plot(wvl1[mask_conv],flux_fit_EW,'r',label=r'Fitted line (g+m+r)'+', vsini='\
                +str(round(popt_i[3],2))+'km/s, vmac='+str(round(popt_i[-1],2))+\
                r'; $\lambda_0=$'+str(round(popt_i[1],3)))
        
    
        plt.title(r'EW$_{\rm i, nc}=$'+str(round(EW_orig,2))+r'$m\AA$, EW$_{\rm i,c}=$'+str(round(EW_convolved,2))+\
                  r'$m\AA$, EW$_{\rm i, gmr}=$'+str(round(EW_fitted,2))+r'$m\AA$')      
        plt.xlabel(r'wavelength ($\AA$)')
        plt.legend()
        plt.show()
        #gaussian_convolved(x,A,lam0,sigma,vsini,beta,vmac,**kwargs)
    if observations: 
        
        A = spec('HD151018')
        wave,flux,b = A.waveflux()
        
        
        line_c,vsin,vmac = 5410.8,67,65
        if fited_directly: 
            fig,axs = plt.subplots(3)
            fig.suptitle("HD 151018: Holgado's values: vsin_GH="+str(vsin)+' km/s; vmac_GH='+str(vmac)+' km/s')
            
            mask = (wave>line_c-4) & (wave<line_c+4)
                    
            pop = curve_fit(gaussian_convolved,wave[mask],flux[mask],\
            p0=[-0.4,line_c,0.1,vsin,1.5,vmac])[0]
                
            flux_fit = gaussian_convolved(wave[mask],*pop)
            flux_fit_GH = gaussian_convolved(wave[mask],-1,pop[1],0.29,vsin,1.5,vmac)
            
            axs[0].plot(wave[mask],flux[mask],label='Observational')
            axs[0].plot(wave[mask],flux_fit_GH,label="Holgados's values")
            axs[0].plot(wave[mask],flux_fit,label='Fitted: vsin='+str(round(pop[-3],2))+' km/s; vmac='+str(round(pop[-1],2))+' km/s')
            axs[0].legend()

            line_c = 5015
    
            mask = (wave>line_c-4) & (wave<line_c+4)
                    
            pop = curve_fit(gaussian_convolved,wave[mask],flux[mask],\
            p0=[-0.4,line_c,0.1,vsin,1.5,vmac])[0]
                
            flux_fit = gaussian_convolved(wave[mask],*pop)
            flux_fit_GH = gaussian_convolved(wave[mask],-0.43,pop[1],pop[2],vsin,1.5,vmac)
            
            axs[1].plot(wave[mask],flux[mask],label='Observational')
            axs[1].plot(wave[mask],flux_fit_GH,label="Gonzalo's values: vsin_GH="+str(vsin)+' km/s; vmac_GH='+str(vmac)+' km/s')
            axs[1].plot(wave[mask],flux_fit,label='Fitted: vsin='+str(round(pop[-3],2))+' km/s; vmac='+str(round(pop[-1],2))+' km/s')
            axs[1].legend()
            
            line_c = 5811.2
    
            mask = (wave>line_c-4) & (wave<line_c+4)
                    
            pop = curve_fit(gaussian_convolved,wave[mask],flux[mask],\
            p0=[-0.4,line_c,0.1,vsin,1.5,vmac])[0]
            flux_fit = gaussian_convolved(wave[mask],*pop)
            flux_fit_GH = gaussian_convolved(wave[mask],-0.33,pop[1],pop[2],vsin,1.5,vmac)
            
            axs[2].plot(wave[mask],flux[mask],label='Observational')
            axs[2].plot(wave[mask],flux_fit_GH,label="Holgados's values")
            axs[2].plot(wave[mask],flux_fit,label='Fitted: vsin='+str(round(pop[-3],2))+' km/s; vmac='+str(round(pop[-1],2))+' km/s')
            axs[2].legend()
            fig.supxlabel(r'Wavelength [$\AA$]')
            
        A.fitline(line_c,func='gc', p0=[-1.0,5410.8,0.25,vsin,1.5,vmac],plot=True)
        A.fitline(5015.0,func='gc', p0=[-0.4,5015.0,0.25,vsin,1.5,vmac],plot=True)
        A.fitline(5811.2,func='gc', p0=[-0.4,5811.2,0.25,vsin,1.5,vmac],plot=True)
    if test:
        
        line_c = 5811.2,5809
        if fited_directly:      
            
            line_c = [5811.2,5809,5806.5]
            fig,axs = plt.subplots(2)
            fig.suptitle("HD 151018: Holgado's values: vsin_GH="+str(vsin)+' km/s; vmac_GH='+str(vmac)+' km/s')
            
            mask = (wave>line_c[0]-6) & (wave<line_c[0]+4)
            pop = curve_fit(gaussian_convolved_general,np.array(wave[mask]),np.array(flux[mask]),\
            p0=[-0.85,-0.25,-0.1,*line_c,0.1,vsin,1.5,vmac])[0]
            #p0=[[-0.85,-0.25],line_c,0.1,vsin,1.5,vmac])[0]
            #ej_res = gaussian_convolved_general(wave[mask],[-0.85,-0.25],line_c,0.1,vsin,1.5,vmac)
            ej_res = gaussian_convolved_general(wave[mask],*pop)
            ej_simple = gaussian_convolved(wave[mask],pop[0],pop[3],*pop[-4:])
            
            axs[0].plot(wave[mask],flux[mask])
            axs[0].plot(wave[mask],ej_res,label='Fitted: vsin='+str(round(pop[-3],2))+' km/s; vmac='+str(round(pop[-1],2))+' km/s')
            axs[0].plot(wave[mask],ej_simple,label='Single line EW='+str(round(sum(ej_simple-1)*np.mean(np.diff(wave[mask]))*-1e3,2))+r'm$\AA$')
            
            axs[0].legend()
            
            line_c = 5410.8,5404.3
            mask = (wave>line_c[0]-9) & (wave<line_c[0]+7)
                    
            pop = curve_fit(gaussian_convolved_general,np.array(wave[mask]),np.array(flux[mask]),\
            p0=[-0.85,-0.25,*line_c,0.1,vsin,1.5,vmac])[0]
            #p0=[[-0.85,-0.25],line_c,0.1,vsin,1.5,vmac])[0]
            #ej_res = gaussian_convolved_general(wave[mask],[-0.85,-0.25],line_c,0.1,vsin,1.5,vmac)
            ej_res = gaussian_convolved_general(wave[mask],*pop)
            ej_simple = gaussian_convolved(wave[mask],pop[0],pop[2],*pop[-4:])
            
            axs[1].plot(wave[mask],flux[mask])
            axs[1].plot(wave[mask],ej_res,label='Fitted: vsin='+str(round(pop[-3],2))+' km/s; vmac='+str(round(pop[-1],2))+' km/s')
            axs[1].plot(wave[mask],ej_simple,label='Single line EW='+str(round(sum(ej_simple-1)*np.mean(np.diff(wave[mask]))*-1e3,2))+r'm$\AA$')  
            axs[1].legend()
        mask = (wave>line_c[0]-9) & (wave<line_c[0]+7)
        result = A.fitline(line_c,p0=[-0.85,-0.2,*line_c,0.1,vsin,1.5,vmac],width=10,func='gc',plot=True)
        flux_line = gaussian_convolved_general(wave[mask], result['A'][0],result['lam0'][0],\
                    result['sigma'],result['vsini'],result['beta'],result['vmac'])
        plt.plot(wave[mask],flux_line,lw=1,c='orange',label='EW='+str(round(sum((flux_line-1)*np.mean(np.diff(wave[mask]))),2)*-1e3))
        plt.legend()
        
        line_c = 5410.8,5404.3
        
        mask = (wave>line_c[0]-9) & (wave<line_c[0]+7)
        result = A.fitline(line_c,p0=[-0.85,-0.2,*line_c,0.1,vsin,1.5,vmac],width=10,func='gc',plot=True)
        flux_line = gaussian_convolved_general(wave[mask], result['A'][0],result['lam0'][0],\
                    result['sigma'],result['vsini'],result['beta'],result['vmac'])
        plt.plot(wave[mask],flux_line,lw=1,c='orange',label='EW='+str(round(sum((flux_line-1)*np.mean(np.diff(wave[mask]))),2)*-1e3))
        plt.legend()