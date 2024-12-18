#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: C. Martinez-Sebastian

GOAL : Abundance measurement. 
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata
from matplotlib.ticker import AutoMinorLocator
import pandas as pd
import pickle
import re
plt.rcParams["mathtext.fontset"] = "cm"
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

sys.path.append('/home/cms/Desktop/PhD/0_PROJECT/pyIACOB-2')
import class_FASTWIND
import spec as sp 
from astroquery.simbad import Simbad
Simbad.add_votable_fields('rv_value')
import FASTWIND_plotting as fw
import spectrum_generation_function as sgf
from db import dir
from tools import romanToInt
trans, cw = fw.considered_transitions()


class abundance_determination():
    def __init__(self,star,Teff,logg):
        '''
        Class for abundance determination

        Parameters
        ----------
        star : str,
            Name of the star to analyze.
            
        Teff : float,
            Effective temperature of the star (in kK)
            
        logg : float,
            log(g); gravity of the star

        Returns
        -------
        None.

        '''
        self.star = star 
        self.star_name = star.replace('HD','HD ')
        self.central_wave_lines = {}
        self.Teff, self.logg = Teff,logg
        for t in trans.values():
            #print(t)
            if t[0]=='1':
                for tran in t[1:]:
                    try: self.central_wave_lines[float(cw[tran])] = tran
                    except: pass

        self.central_wave_lines[4554.00 ] = 'SiIII' # As we are constantly seeing Si III close
        self.central_wave_lines[4631.24] = 'SiIV' # As we are constantly seeing Si IV close
        self.central_wave_lines[6378.34] = 'O3' # As we are constantly seeing Si III close # Transition name
        self.central_wave_lines[6375] = 'unk'
        self.central_wave_lines[4515.811] = 'C3'
        self.central_wave_lines[4621.249] = 'O2'
        #central_wave_lines[6379.64] = 'unk' # Is (most probable) O3. ()
        # (and confusing with N II 4552) we added this line for clarifying
        self.central_wave_lines_sorted = np.array(sorted(self.central_wave_lines))
        return None

    def expected_ion(self):
        '''
        Function as a first approximation to knwo, as a function of temperature, which
        C, N, O ions we would expect
        
        Parameters
        ----------
        self,
            
        Returns
        -------
        CII,CIV,NII,NIII,NIV,NV ,OII,OIII,OIV,OV : bool
            If we expect the ion is True
        '''
        CII,CIII,CIV = False,False,False
        NII,NIII,NIV,NV = False,False,False,False
        OII,OIII,OIV,OV = False,False,False,False
        if self.Teff<=29: CII = True
        if self.Teff<=35: NII = True
        if self.Teff>=35: CIV, NIV = True,True
        if self.Teff<36: OII = True
        if self.Teff>43: NV = True
        if self.Teff>45: OV = True
        if self.Teff<50: 
            CIII,NIII = True,True
            if self.Teff>280: OIII = True
            if self.Teff>300: OIV = True
        return CII,CIII,CIV,NII,NIII,NIV,NV,OII,OIII,OIV,OV

    def micro_func(self):
        '''
        Function to give a proxy of the expected microturbulence according to Armin's 
        results
    
        Parameters
        ----------
        None
    
        Returns
        -------
        None
        '''
        self.Lspec = 4 * np.log10( self.Teff * 1000) - self.logg - 10.61
        print(self.Lspec)
        self.mic_g = -11.54 * self.logg + 49.6
        if self.Lspec < 3.5: 
           self.mic_L = 4
        else:
            self.mic_L = 29.3 * self.Lspec - 97.7
        return None
    
    def line_measurement(self,lines,save=True,fitfunc='g',path = dir().modeldir,
                         vsini=1e-3,vmac=1e-3,elem='N',save_csv=True,extra='',txt=False,file=None):
        '''
        Function to study the different lines, if expected. 
        Parameters
        ----------
            
        lines : dict or pd.data,
            Dictionary containing lines to be studied in the star. They could be read
            from 'Ostars_N_isolated_lines.csv' or equivalent. The structure is: 
            lines = {"LINES":[...],star_name:[line_activation_values]}
            
        save : bool, optional   
            If True, it saves the measurement figures. Default is True. 
            
        fitfunc : str, optional
            Determines the function spec is using for line fitting. Default is 'g' (Gaussian)
            Accepts also: 'gc' (Gaussian covolved with radial-tangencial profile+rotation),
            'l' (lorentzian), 'v' (Voigt profile), 'r' (rot)...
            
        path : str, optional
            Path to save images (followed by "IMAGES/*elem*_DETERMINATION/") and csv
            files (followed by "star_analysis"). Default is dir().modeldir
     
        vsini : float, optional
            Rotational velocity. First approach to use (if using covolved profile). 
            Default is 1e-3
            
        vmac : float, optional
            Macroturbulence velocity. First approa
            
        elem : str, optional
            Chemical symbol of the element which lines are to be stuied. 
            Determines folder to save images + csv extention in star_abalysis.
            Default is "N"
            
        save_csv : bool, optional
            If True, the results (dict_star) are saved into a csv. Default is True.
        
        extra: str, optional
            Extra string chain to add in the newly created files, when necessary. 
            Default is ''.
            
        txt : boolean, optional
            If True, it assumes spectrum from a two-columns file with wavelenght and flux
            with no header in the file. Default is False.
            
        file : str, optional
            Name of the spectra to be used (including path). If none is given, assumes
            star name and looks by default in spectra direction. Default is None.
            
        Returns
        -------
        dict_star : dict, 
            Dictionary containing the different lines, activation value (1 if EW is reliable, 0 if not), EW, and comments. 
    
        '''
        dict_star = {}
        lambda_0,Comments,EW_obs_list,Activation = {},{},{},{}
        if file is None: file = self.star
        ##################### CHARGING OBSERVATIONAL SPECTRUM #####################
        q = Simbad.query_object(self.star)
        if self.star == 'HD96946': q['RV_VALUE'][0] = 0.
        #elif star == 'HD48279': q['RV_VALUE'][0] = -400.
        # Correcting from velocity variation 
        try: 
            if type(q['RV_VALUE'][0])==float:
                RV_compute = sp.spec(file,SNR='best',rv0=q['RV_VALUE'][0],txt=txt)
                print('RV value:',q['RV_VALUE'][0])
            else: RV_compute = sp.spec(file,SNR='best',txt=txt)
        except:  RV_compute = sp.spec(file,SNR='best',txt=txt)
        
        print(q['RV_VALUE'][0])
        fit_RV = RV_compute.fitline(5592.37,plot=True,p0=[-0.08,5592.37,0.5])
        #plt.plot(np.linspace(5586,5600,1000), sp.f_gaussian1(np.linspace(5586,5600,1000),-0.08,5592.37,0.5))    
        if fit_RV['sol'] == 0: fit_RV = RV_compute.fitline(5592.37,plot=True,p0=[-0.04,5592.37,1.5]);print('Fitting broad line 5592')
        if fit_RV['sol'] == 0: fit_RV = RV_compute.fitline(4552.622,plot=True,p0=[-0.15,4552.622,0.5]);print('Fitting line 4552')
        if fit_RV['sol'] == 0: fit_RV = RV_compute.fitline(4603.73 ,plot=True,p0=[-0.15,4603.73,0.5]);print('Fitting line 4603')
        if fit_RV['sol'] == 0: 
            print('No initial fit')
            return None
        try: 
            if type(q['RV_VALUE'][0])==float and fit_RV['sol']==1:
                observational_spectrum = sp.spec(file,SNR='best',rv0 = fit_RV['RV_kms']+q['RV_VALUE'][0],txt=txt)
            elif fit_RV['sol']==1:
                observational_spectrum = sp.spec(file,SNR='best',rv0 = fit_RV['RV_kms'],txt=txt)
        except: observational_spectrum = sp.spec(file,SNR='best',rv0 = fit_RV['RV_kms'],txt=txt)
        
        FWHM_ref = fit_RV['FWHM']
        print('REFERENCE VALUE FWHM: %.2f' % (FWHM_ref))
        if txt: observational_spectrum.txtwaveflux()
        else: observational_spectrum.waveflux()
        ##################### CHARGING OBSERVATIONAL SPECTRUM #####################
        
        
        if self.star not in lines:
            print(self.star+' not in lines')
            lines[self.star] = np.ones(len(lines['LINES']))
        # Reading the lines to consider
        if save_csv:
            lines.to_csv('Ostars_'+elem+'_isolated_lines.csv', index=False)
                
        star_lines = dict(zip(lines["LINES"], lines[self.star]))
        
        # Generating a multi-page figure to see the fitting of the different lines
        if save and save_csv:
            if elem == 'N':
                pp = PdfPages(path+'/IMAGES/NITROGEN_DETERMINATION/'+self.star+'_'+fitfunc+extra+"_multi.pdf")
            if elem == 'C':
                pp = PdfPages(path+'/IMAGES/CARBON_DETERMINATION/'+self.star+'_'+fitfunc+extra+"_multi.pdf")
            if elem == 'O':
                pp = PdfPages(path+'/IMAGES/OXYGEN_DETERMINATION/'+self.star+'_'+fitfunc+extra+"_multi.pdf")
        elif save:
            if elem == 'N':
                pp = PdfPages(path+'/IMAGES/NITROGEN_DETERMINATION/'+self.star+'_'+fitfunc+extra+"_multi_add.pdf")
            if elem == 'C':
                pp = PdfPages(path+'/IMAGES/CARBON_DETERMINATION/'+self.star+'_'+fitfunc+extra+"_multi_add.pdf")
            if elem == 'O':
                pp = PdfPages(path+'/IMAGES/OXYGEN_DETERMINATION/'+self.star+'_'+fitfunc+extra+"_multi.pdf")
        # See if the different ions are expected. 
        C2,C3,C4,N2,N3,N4,N5,O2,O3,O4,O5 = self.expected_ion()
        HE, H1 = True, True
        Si = True
        un = True
        for lin in lines['LINES']: 
            print('')
            print(lin)
            elem = lin[0]
            central_wave = float(cw[trans[lin][1]])
            norm=True 
            ######### LOONKING FOR THE DIFFERENT LINES AROUND THE STUDIED #########
            MIN,MAX = central_wave,central_wave
            
            # Lines around the studied
            lines_to_consider1 = self.central_wave_lines_sorted[(self.central_wave_lines_sorted>MIN-2*FWHM_ref) & (self.central_wave_lines_sorted<MAX+2*FWHM_ref)]
            lines_to_consider1
            # Studying if those lines are expected.
            ltc1,ltc2 = [],[]
            for l in lines_to_consider1:
                if self.central_wave_lines[l][:2] != '23':
                    if locals()[self.central_wave_lines[l][:2]]: ltc1.append(l)
                else:
                    if N2: ltc1.append(l)
            if central_wave not in ltc1: ltc1.append(central_wave)
            lines_to_consider2 = self.central_wave_lines_sorted[(self.central_wave_lines_sorted>min(ltc1)-2*FWHM_ref)\
                                                           & (self.central_wave_lines_sorted<max(ltc1)+2*FWHM_ref)]    
            # If the current lines considered are less than increasing the area, repeat the study
            if len(lines_to_consider1)!=len(lines_to_consider2):
                for l in lines_to_consider2:
                    if self.central_wave_lines[l][:2] != '23':
                        print(locals())
                        if locals()[self.central_wave_lines[l][:2]]: ltc2.append(l)
                    else:
                        if N2: ltc2.append(l)
                if central_wave not in ltc2: 
                    ltc2.append(central_wave)
                while len(ltc1) != len(ltc2):
                    ltc1 = ltc2
                    ltc2 = []
                    lines_to_consider2 = self.central_wave_lines_sorted[(self.central_wave_lines_sorted>min(ltc1)-2*FWHM_ref)\
                                                                   & (self.central_wave_lines_sorted<max(ltc1)+2*FWHM_ref)]
                    for l in lines_to_consider2:
                        if self.central_wave_lines[l][:2] != '23':
                            try:
                                if locals()[self.central_wave_lines[l][:2]]: ltc2.append(l)
                            except: print('problems in',l)
                        else:
                            if N2: ltc2.append(l)  
            # Adding the line to study if it has been removed
            if central_wave in ltc1: 
                ltc1.remove(central_wave)
            ltc1.insert(0, central_wave)
            ######### LOONKING FOR THE DIFFERENT LINES AROUND THE STUDIED #########
            
    
            ### DETERMINING INITIAL CONDITIONS AND LIMITS FOR THE LINES TO CONSIDER ###
            p,bound_min,bound_max =  [],[],[]
            if fitfunc == 'g':
                for c_w in ltc1:
                    mean_flux = np.mean(observational_spectrum.flux[(observational_spectrum.wave>c_w-FWHM_ref/2) &\
                                            (observational_spectrum.wave<c_w+FWHM_ref/2)])
                    if mean_flux<1:
                        p.append(min(observational_spectrum.flux[(observational_spectrum.wave>c_w-FWHM_ref/2) &\
                                                (observational_spectrum.wave<c_w+FWHM_ref/2)])-1)
                        bound_min.append(-1),bound_max.append(0)
                    else:
                        p.append(max(observational_spectrum.flux[(observational_spectrum.wave>c_w-FWHM_ref/2) &\
                                                (observational_spectrum.wave<c_w+FWHM_ref/2)])-1)
                        bound_min.append(0),bound_max.append(1)
                    p.append(c_w),bound_min.append((1-1/3e4)*c_w),bound_max.append((1+1/3e4)*c_w)
                    p.append(FWHM_ref/2.355),bound_min.append(0.85*FWHM_ref/2.355),bound_max.append(1.15*FWHM_ref/2.355)
    
            elif fitfunc == 'gc':
                p_A,p_lam = [],[]
                bound_min.extend((-1*np.ones(len(ltc1)))),bound_max.extend((1*np.ones(len(ltc1))))
                
                for c_w in ltc1:
                    bound_min.append(c_w-(FWHM_ref/2.355)),bound_max.append(c_w+(FWHM_ref/2.355))
                    EW_proxy = sum(observational_spectrum.flux[(observational_spectrum.wave>c_w-FWHM_ref) &\
                                            (observational_spectrum.wave<c_w+FWHM_ref)]-1)*np.mean(np.diff(observational_spectrum.wave))
                    A_proxy = np.sqrt(2*np.pi)*EW_proxy/np.pi/(FWHM_ref/2.355)
                    print('A proxy is ',A_proxy)
                    p_A.append(A_proxy)
                    p_lam.append(c_w)
                bound_min.extend([0,0,0,0]),bound_max.extend([3,200,1.5,200])
                p = *p_A,*p_lam,0.1,vsini,1.5,vmac
                # Fitting the gaussians 
            
            bounds = (bound_min,bound_max)
            print('Initial values are',p,'and bounds',bounds)    
            
            # Fittinf the line
            fit = observational_spectrum.fitline(ltc1,line=central_wave,plot=True,\
            lim=bounds,p0=p,normalization=norm,sigma_tres=2.,func=fitfunc)
                
            # If we are satisfied with the solution, we save it. 
            try: fit['lines'] = np.array(fit['lines']).astype('float32')
            except: pass
        
            if lin in EW_obs_list.keys() and fit['sol']==0:
                pass
            elif fit['sol']==0 or round(central_wave,2) not in np.around(fit['lines'],2) or fit['EW']<15 or star_lines[lin] == 0:
                
                Activation[lin] = 0
                lambda_0[lin] = central_wave
                #EW_obs_list.append(fit['EW'])
                if lin not in EW_obs_list.keys() or EW_obs_list[lin] == None:
                    EW_obs_list[lin] = fit['EW']
                Comments[lin] = 'Non detected'
                star_lines[lin] = 0
                try:
                    plt.title(lin+', '+str(central_wave)+' not considered; EW='+str(fit['EW']))
                except: 
                    plt.title(lin+', '+str(central_wave)+' not considered')
            elif  lin != 'NIII4511' and lin != 'NIII4515' and fit['FWHM']>5:
                Activation[lin] = 0
                lambda_0[lin] = central_wave
                #EW_obs_list.append(fit['EW'])
                if lin not in EW_obs_list.keys() or EW_obs_list[lin] == None:
                    EW_obs_list[lin] = fit['EW']
                Comments[lin] = 'Non detected'
                star_lines[lin] = 0
                try:
                    plt.title(lin+', '+str(central_wave)+' not considered; EW='+str(fit['EW']))
                except: 
                    plt.title(lin+', '+str(central_wave)+' not considered')
            else:
                Activation[lin] = 1
                lambda_0[lin] = central_wave
                # EW_obs_list.append(fit['EW_dict'][fit['line']])
                if lin not in EW_obs_list.keys() or EW_obs_list[lin] == None:
                    EW_obs_list[lin] = fit['EW_dict'][fit['line']]
                Comments[lin] = ''
                try:
                    plt.title(lin+', '+str(central_wave)+'; EW='+str(fit['EW_dict'][fit['line']])\
                    +r' m$\AA$; vsini='+str(fit['vsini'])+' km/s; vmac='+str(fit['vmac'])+' km/s')
                except:
                    plt.title(lin+', '+str(central_wave)+'; EW='+str(fit['EW']))
                    
                EW_obs = fit['EW']
                #if lin=='NIII4518':
                #    print('fit is:',fit['EW_dict'].keys(),type(fit['EW_dict'].keys()))
                #    break
                for critical_line in '4510.965','4514.817','4514.854','4518.143','4523.561':
                    print('WARNING: critical lines: ',critical_line,round(float(critical_line),3) in fit['EW_dict'].keys(),fit['EW_dict'].keys(),type(fit['EW_dict'].keys()))
                    if round(float(critical_line),3) in fit['EW_dict'].keys():
                        if critical_line == '4510.965': line2 = 'NIII4511'
                        elif critical_line == '4514.817' or critical_line == '4514.854': line2 = 'NIII4515'
                        elif critical_line == '4518.143': line2 = 'NIII4518'
                        elif critical_line == '4523.561': line2 = 'NIII4523'
                        print('WE ARE IN CRITICAL CASE: ',line2)
                        if line2 not in EW_obs_list.keys() or EW_obs_list[line2] == None:
                            EW_obs_list[line2] = fit['EW_dict'][float(critical_line)]
                            lambda_0[line2] = central_wave
                            if EW_obs_list[line2] >= 15:
                                Activation[line2] = 1
                                Comments[line2] = ''
                            else:
                                Activation[line2] = 0
                                Comments[line2] = 'Non detected'     
            # Adding the desired information to the fitting figures. 
            fig_nums = plt.get_fignums()
            fig = plt.figure(fig_nums[-1])
            fig_nums = plt.get_fignums()    
            for c in ltc1:  #### QUE NOS DIGA CUÁLES SE USAN Y CUÁLES NO 
                txt = plt.annotate(self.central_wave_lines[c],xy=(c,1.))
                txt.set_rotation(45)
                plt.axvline(c,c='r',ls='--')
                #print('FOR ',lin,', VERTICAL LINES ARE',c)
            try: # Plotting the considered lines
                if fit['EW']>=15:
                    if fitfunc == 'g':
                        for f in range(len(fit['lines'])//3):
                            plt.axvline(fit['lines'][int(3*f+1)],c='b',ls='--')
                    elif fitfunc == 'gc':
                        for f in range(len(fit['lines'])//2-2):
                            plt.axvline(fit['lines'][len(fit['lines'])//2-2+f],c='b',ls='--')
        
            except: print('lines do not exist for',lin)
            
            if save:
                fig.savefig(pp, format='pdf')
        if save:
            pp.close()   
                    
            plt.close('all')
        print(len(EW_obs_list.keys()),len(Activation.values()),len(lambda_0.values()),\
              len(EW_obs_list.values()),len(Comments.values()))
        # Creating and saving a dictionary for the fitting results of the studied star
        dict_star['line_name'],dict_star['Activation'],dict_star['lambda_0'],dict_star['EW_obs'],dict_star['Comments'] = \
        EW_obs_list.keys(),Activation.values(),lambda_0.values(),EW_obs_list.values(),Comments.values()
        if save_csv:
            df = pd.DataFrame(dict_star)
            df.to_csv(path+'/star_analysis/'+self.star+'_'+elem+'_'+fitfunc+extra+'.csv',index=False) # We added an N to diferenciate: it only has Nitrogen lines. 
        print('DIRECTORY: ',dict_star)
        return dict_star

    @staticmethod
    def D_EW(line,model,EW_obs,dEW = False,step_number_abundance = 420,step_number_microturbulence = 30,\
             plot2D = False, plot3D = False):
        '''
        Function to interpolate the abundance and microturblence and study the abubndance 
        in two dimessions. It uses the EW from FASTWIND model directly. 
        
        It uses chi2 = ((EW_obs - EW_mod) / EW_obs)**2 * dEW
    
        Parameters
        ----------
        line : str,
            Name of the line to be studied.
            
        model : class Model,
            Element model class.
            
        EW_obs : float,
            Value of the measured EW for the observational spectrum.
            
        dEW : float, optional
            Uncertainty of the measured EW (in mAA). If False, it is taken as 10% of 
            the EW_obs. Default is False.
            
        step_number_abundance : int, optional
            Number of steps for the abundance to be interpolated to from the minimum
            to the maximum form the models. The default is 105.
            
        step_number_microturbulence : int, optional
            Number of steps for the microturbulence to be interpolated to from the 
            minimum to the maximum form the models. The default is 30.
            
        plot2D : bbol, optional
            If true, it generates a density map with the results (abundace, microturbulence
            and chi2_EW). The default is False.
            
        plot3D : bool, optional
            If true, it generates a surface with the results (abundace, microturbulence
            and chi2_EW). The default is False.
    
        Returns
        -------
        ab,mic,chi2_EW
        
        ab_linearized : np.array
            Interpolated abundance.
            
        mic_linearized : np.array
            Interpolated microturbulence.
            
        chi2_linearized : np.array
            Resulting chi2 for mapping the EW obs respect to FW.
        '''
        # Dictionary of the line to be generated
        model_line = {'abundance' : [], 'microturbulence' : [], 'EW' : []}
        ################ Generating the original model dictionary ################
        for ab in model.lines[line].abundances.keys(): # Computing pseudo-xi
            for mic in model.lines[line].abundances[ab].microturbulences.keys():
                try: 
                    model_line['EW'].append(model.get_profile(line,ab,mic).EW)#*1e6#/EW_obs)) # We multiply by 1e6 due to the error when saving EW
                    model_line['abundance'].append(float(ab)+12) # epsilon_N
                    model_line['microturbulence'].append(float(mic[-3:]))
                except: pass#print(ab,mic)
        if not dEW: dEW = 0.1 * EW_obs
        ab,mic,ew = np.array(model_line['abundance']),np.array(model_line['microturbulence']),np.array(model_line['EW'])
        ew[ew==None] = 1e-6
        sim_v = max(ew)
        ########## Interpolation ##########
        ab_lin = np.linspace(min(np.array(model_line['abundance'])),max(np.array(model_line['abundance'])),step_number_abundance)
        mic_lin = np.linspace(min(np.array(model_line['microturbulence'])),max(np.array(model_line['microturbulence'])),step_number_microturbulence)
        
        AB,MIC = np.meshgrid(ab_lin,mic_lin)
        EW = griddata((ab,mic), ew, (AB,MIC), method='cubic')   
        chi2 = ((EW_obs-EW)/EW_obs)**2*dEW
        
        ab_linearized = AB.reshape(len(ab_lin)*len(mic_lin))
        mic_linearized = MIC.reshape(len(ab_lin)*len(mic_lin))
        chi2_linearized = chi2.reshape(len(ab_lin)*len(mic_lin))
        ########## Interpolation ##########
        
        if plot3D:
            ########## 3D plot ##########
            fig = plt.figure()
            ax1 = fig.add_subplot(121, projection = '3d')
            ax2 = fig.add_subplot(122, projection = '3d', sharex = ax1, sharey = ax1, sharez = ax1)
            chi2 = (1-ew)**2
            contourf_ = ax1.plot_trisurf(ab,mic,chi2,vmax=sim_v,vmin=0,cmap = cm.afmhot)
            ax2.plot_surface(AB,MIC,chi2,vmax=sim_v,vmin=0,cmap = cm.afmhot)
            
            ax1.set_xlabel(r'$\epsilon$'),ax1.set_ylabel(r'$\xi_t$ [km/s]'),ax1.set_zlabel(r'$\frac{EW_{mod}-EW_{obs}}{EW_{obs}}$')
            ax2.set_xlabel(r'$\epsilon$'),ax2.set_ylabel(r'$\xi_t$ [km/s]'),ax2.set_zlabel(r'$\frac{EW_{mod}-EW_{obs}}{EW_{obs}}$')
            fig.subplots_adjust(right=0.8)
            cbar_ax = fig.add_axes([0.9, 0.15, 0.03, 0.7])
            fig.colorbar(contourf_, cax=cbar_ax)
            fig.suptitle(line)
            ########## 3D plot ##########
        
        if plot2D:
            ########## Contour plot ##########
            chi2 = (1-ew)**2
            ab_sorted,mic_sorted,chi2_sorted = zip(*sorted(zip(ab,mic,chi2)))
            n,m = len(np.unique(ab)),len(np.unique(mic))
            chi2_2D = np.array(chi2_sorted).reshape(n,m)
            ab_2D,mic_2D = np.array(ab_sorted).reshape(n,m), np.array(mic_sorted).reshape(n,m)
            
            fig,(ax1,ax2) = plt.subplots(2)
            contourf_ = ax1.pcolormesh(ab_2D,mic_2D,chi2_2D,vmax=sim_v,vmin=0,cmap = cm.afmhot)
            ax2.pcolormesh(ab_2D,mic_2D,chi2_2D,vmax=sim_v,vmin=0,cmap = cm.afmhot)
            
            ax1.set_ylabel(r'$\xi_t$ [km/s]'),ax2.set_xlabel(r'$\epsilon$'),ax2.set_ylabel(r'$\xi_t$ [km/s]')
            fig.suptitle(line)
            
            fig.subplots_adjust(right=0.8)
            cbar_ax = fig.add_axes([0.9, 0.15, 0.03, 0.7])
            fig.colorbar(contourf_, cax=cbar_ax)
            ########## Contour plot ##########
            plt.show()
        return ab_linearized,mic_linearized,chi2_linearized

    def EW_study(self, line_weights = False):
        '''
        Pure EW study. 
        Here, resulting abundance is weighted by the EW and EW uncertainty.
        
        Parameters
        ----------
        line_weights : dict, optinal
            Dictionary (same key-names as elem_list) with the weight of each line. 
            If False, each line if weighted as 1. Deffault is False
    
        Returns
        -------
        None
        '''
        elem = list(self.elem_list.keys())[0].replace('I','').replace('V','')
        self.mic_FWgrid = pd.read_csv(dir().modeldir+'grid_components/microturbulences.dat', header=None,sep='\s+')[0]
        self.ab_FWgrid = pd.read_csv(dir().modeldir+'grid_components/abundances_'+elem+'.dat', header=None,sep='\s+')[1]
        i,j = 0,0 # weight of TOTAL lines, number of ions (with elements)
        elem_EW_chi2 = 0
        ab,mic,EW_chi2 = {},{},{}
        N_ion = {}
        for ion in self.elem_list.keys():
            if len(self.elem_list[ion]) != 0:
                j += 1
                k = 0 # weight of ION lines
                ab[ion],mic[ion],EW_chi2[ion] = 0,0,0
                for lin in self.elem_list[ion].keys():
                    #if line_weights: w = 1/line_weights[lin] # number of lines with this transition
                    #else: w = 1
                    w = 1
                    ab[ion] = self.elem_list[ion][lin]['abundances']
                    mic[ion] = self.elem_list[ion][lin]['microturbulences']
                    chi2 = self.elem_list[ion][lin]['chi2_EW']
                    max_val = sorted(chi2)[len(chi2)//10+1] # Chosing only 10% of best fitting values. Worse, settle to 1. 
                    last_val = sorted(chi2)[len(chi2)//10]/max_val
                    chi2 /= max_val # Normalizing to the best 10 %
                    chi2[chi2>last_val] = 1
                    #chi2 = (chi2+ w -1) / w
                    chi2 /= w
                    EW_chi2[ion] += chi2
                    i += 1/w
                    k += 1/w
                    self.elem_mic,self.elem_ab  = mic[ion],ab[ion]
                N_ion[ion] = round(k,2)
                ab[ion],mic[ion],EW_chi2[ion] =\
                ab[ion],mic[ion],EW_chi2[ion]
                elem_EW_chi2 += EW_chi2[ion]
                EW_chi2[ion] = EW_chi2[ion]/k
        if i==0:
            print('WARNING: Element %(elem)s has an empty list.' % {"elem": elem})
            return None, None, None
        self.elem_EW_chi2 = elem_EW_chi2/i
        print('Maximum value is',max(self.elem_EW_chi2),'(should be 1)')
        n,m = len(np.unique(self.elem_mic)),len(np.unique(self.elem_ab))
        
        # Creating summary figure for each ion (and the general one, summing all of them)
        fig = plt.figure()
        ax = {}
        if j%2 == 0:
            for i in range(j):
                ax[i] = fig.add_subplot(j//2+1,2,i+1)
            ax[i+1] = fig.add_subplot(j//2+1,1,j//2+1)
        else:
            for i in range(j):
              ax[i] = fig.add_subplot(j//2+1,2,i+1)
            ax[i+1] = fig.add_subplot(j//2+1,2,i+2)  
        i = 0
        for ion in self.elem_list.keys():
            for mic_grid_tick in list(self.mic_FWgrid):
                ax[i].axhline(float(mic_grid_tick.replace('VT','')),lw=1,c='grey',alpha=0.8,zorder=100,ls=':')
            for ab_grid_tick in list(self.ab_FWgrid):
                ax[i].axvline(ab_grid_tick+12,lw=1,c='grey',alpha=0.8,zorder=100,ls=':')
            if len(self.elem_list[ion]) != 0:
                fw.plot_style(ax[i],size=12)
                max_level = max(EW_chi2[ion][EW_chi2[ion]<1])
                levels = np.linspace(min(EW_chi2[ion]), max_level, 21)
                ax[i].contourf(ab[ion].reshape(n,m),mic[ion].reshape(n,m),EW_chi2[ion].reshape(n,m),\
                            levels=levels,cmap = cm.YlOrBr_r)
                    
                decim = -int(np.log10(min(EW_chi2[ion])))+1
                if decim<0: decim=2
                
                ax[i].plot(ab[ion][EW_chi2[ion]==min(EW_chi2[ion])],\
                    mic[ion][EW_chi2[ion]==min(EW_chi2[ion])],'rx',\
                    label=r'$\epsilon$='+str(round(ab[ion][EW_chi2[ion]==min(EW_chi2[ion])][0],2))+\
                    r', $\xi_t$='+str(round(mic[ion][EW_chi2[ion]==min(EW_chi2[ion])][0],2))+\
                    r', $\chi^2=$'+str(round(min(EW_chi2[ion]),decim)))
        
                ax[i].plot([], [], ' ', label='$N_{\\rm lines}=$'+str(len(self.elem_list[ion])))
                ax[i].set_title(ion)
                ax[i].legend()
                ax[i].set_xlim(min(self.ab_FWgrid)+12,max(self.ab_FWgrid)+12); ax[i].set_ylim(float(min(self.mic_FWgrid).replace('VT','')),float(max(self.mic_FWgrid).replace('VT','')))
    
                i += 1
        fw.plot_style(ax[i],size=12)
        max_level= max(self.elem_EW_chi2[self.elem_EW_chi2<1])
        levels = np.linspace(min(self.elem_EW_chi2), max_level, 21)        
        contourf_el = ax[i].contourf(self.elem_ab.reshape(n,m),self.elem_mic.reshape(n,m),self.elem_EW_chi2.reshape(n,m),\
                    levels=levels,cmap = cm.YlOrBr_r)  
        decim = -int(np.log10(min(self.elem_EW_chi2)))+1
        if decim<0: decim=1
    
        for mic_grid_tick in list(self.mic_FWgrid):
            ax[i].axhline(float(mic_grid_tick.replace('VT','')),lw=1,c='grey',alpha=0.8,zorder=100,ls=':')
        for ab_grid_tick in list(self.ab_FWgrid):
            ax[i].axvline(ab_grid_tick+12,lw=1,c='grey',alpha=0.8,zorder=100,ls=':')
        ax[i].set_xlim(min(self.ab_FWgrid)+12,max(self.ab_FWgrid)+12); ax[i].set_ylim(float(min(self.mic_FWgrid).replace('VT','')),float(max(self.mic_FWgrid).replace('VT','')))
        print(self.elem_ab[np.argmin(self.elem_EW_chi2)])
        ax[i].plot(self.elem_ab[np.argmin(self.elem_EW_chi2)], self.elem_mic[np.argmin(self.elem_EW_chi2)],'rx',\
             label=r'$\epsilon$='+str(round(self.elem_ab[np.argmin(self.elem_EW_chi2)],2))+\
             r' $\xi_t$='+str(round(self.elem_mic[np.argmin(self.elem_EW_chi2)],2))+r', $\chi^2=$'+str(round(min(self.elem_EW_chi2),\
             decim)))
        ax[i].legend()
        ax[i].set_title(ion[0])
        # Colorbar
        fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.85, 0.15, 0.03, 0.7])
        cbar = fig.colorbar(contourf_el, cax=cbar_ax, ticks=[0,max_level])
        cbar.ax.set_yticklabels(['0', 'best 10%']) 
        cbar.set_label(r'$\chi^2=\left(\sum w\right)^{-1}\sum_{\rm line}\left(\frac{EW_{\rm mod}-EW_{\rm obs}}{EW_{\rm obs}}\right)^2\Delta EW_{\rm obs}w$',loc='center', fontsize=14,labelpad=-20)
        fig.supxlabel(r'$\log(\frac{\rm '+ion[0]+'}{\\rm H})+12$ [dex]', fontsize=16),fig.supylabel(r'$\xi_t$ [km/s]', fontsize=12)
        plt.show()
        return None#elem_mic,elem_ab,elem_EW_chi2

    def abundance_determination_func(self, elem, model, vsini, vmac, SpC='Unknown',
                                observed_EW_dic=False, fix_micro=False, mic_L = False,
                                scatter_minimum=True, step_number_abundance=420,
                                step_number_microturbulence=30, save_im = False,
                                save_im_path=dir().modeldir+'IMAGES/', extra='',
                                FW_model_path=dir().modelgrid+dir().modelgrid.split('/')[-2]+'_pkl/original/'):
        '''
        Final function to stimate micro and abundance from different ways. 
        
        Parameters
        ----------
        elem : str
            Element to compute the abundance. Is used to save the information and 
            compute the ionization state.
        model : str,
            FASTWIND model name (should be in the folder FW_model_path). Should 
            have the format "TXXXgXXXHeXXQXXXXbXX".
        vsini : float
            Projected rotational velocity of the star. Could be obtaines in 
            IACOB-Broad.
        vmac : float
            Macroturbulence velocity of the star. Could be obtaines in IACOB-Broad.
        SpC : str, optional
            Spectral type of the star for titles.
            The default is 'Unknown'.
        observed_EW_dic : pd.DataFrame, optional
            Pandas dataframe with line-name in 'line_name', EW in 'EW_obs', and 
            aconsideredness in 'Acivation' (1 if considered, 0 otherwise).
            If not given, it read it in csv 'dir().modeldir+'star_analysis/'+star+'_'+elem+'_gc.csv'
        fix_micro : float, optional
            If given, it use tis value as micro. 
            The default is False.
        mic_L : bool, optional
            If True, it stimates the micro with luminosity using 'micro_func'.
            The default is False.
        scatter_minimum : bool, optional
            If True, the micro if fixed to the value that minimize the scatter.
            The default is True.
        step_number_abundance : int, optional
            Number of abundance steps to interpolate the grid to. 
            The default is 420.
        step_number_microturbulence : int, optional
            Number of microturbulence steps to interpolate the grid to. 
            The default is 30.
        save_im : bool, optional
            If True, it saves the images. 
            The default is False.
        save_im_path : str, optional
            Path in which save images. 
            The default is dir().modeldir+'IMAGES/'.
        extra : str, optional
            String chain to add at the end of figure saving names. 
            The default is ''.

        Returns
        -------
        {'Micro', 'Weighted abundance', 'Weighted uncertainty'}.
        '''
        self.elem = elem
        self.model, self.vsini, self.vmac, self.SpC = model, vsini, vmac, SpC
        
        self.mic_min_sct_ion,self.mic_min_sct_multiplet = {},{}
        if observed_EW_dic is False:
            observed_EW_dic = pd.read_csv(dir().modeldir+'star_analysis/'+star+'_'+self.elem+'_gc.csv')
        self.elem_list = {} # elem_list[ion][lin] = {'abundances', 'microturbulences','chi2_EW'}.

        #self.abundances,self.uncertainties,self.microturbulences = {},{},{}
        
        ### Charging FW model ###
        with open(FW_model_path+self.model+elem+'_lines.pkl','rb') as inp:
            self.FW_model = pickle.load(inp, encoding='latin1')
        ### Charging FW model ###
        
        ### Weighting the lines (avoiding overconsideration of transitions) ###
        tran_dic,self.line_weight = {},{}
        for l in range(len(observed_EW_dic['line_name'])):
            line = observed_EW_dic['line_name'][l]
            if observed_EW_dic['Activation'][l]==1:
                transitions = trans[line]
                for i in range(1,len(transitions)):
                    if transitions[i][:-2] not in tran_dic.keys():
                        tran_dic[transitions[i][:-2]] = [line]
                    elif line not in tran_dic[transitions[i][:-2]]: 
                        tran_dic[transitions[i][:-2]].append(line)
        for l in range(len(observed_EW_dic['line_name'])):
            line = observed_EW_dic['line_name'][l]
            if observed_EW_dic['Activation'][l]==1:
                transitions = trans[line]
                n = 0
                for i in range(1,len(transitions)):
                    if len(tran_dic[transitions[i][:-2]]) > n: n = len(tran_dic[transitions[i][:-2]])
                #print(line,1/n)
                self.line_weight[line] = 1/n
        ### Weighting the lines (avoiding overconsideration of transitions) ###
        CII,CIII,CIV,NII,NIII,NIV,NV,OII,OIII,OIV,OV = self.expected_ion()
        for l in range(len(observed_EW_dic['line_name'])):
            if observed_EW_dic['Activation'][l]==1 and observed_EW_dic['line_name'][l].startswith('N'):
                lin = observed_EW_dic['line_name'][l]
                
                ion = lin[:re.search(r"\d",lin).start()]
                # Mapping in 2d the EW difference between observations and model FOR A GIVEN LINE 
                ab,mic,chi2_EW = self.D_EW(lin,self.FW_model,observed_EW_dic['EW_obs'][l]) # This will be latter used
                if ion not in self.elem_list.keys():
                    self.elem_list[ion] = {}
                if self.FW_model:
                    # Saving the EW difference for each line (divided by ions)
                    self.elem_list[ion][lin] = {'abundances' : ab, 'microturbulences' : mic, 'chi2_EW' : chi2_EW}    
                        
        ### FIRST APPROACH TO THE ABUNDANCE: directly mapping all lines (EW study and plots)
        self.EW_study(line_weights = self.line_weight)
   
        plt.suptitle(self.star_name+', '+self.model+'\n' '$v\sin(i)$='+str(self.vsini)+r' [km/s], v$_{\rm mac}$='+str(self.vmac)+' [km/s]')
        if save_im: plt.savefig(save_im_path+'Density_Map_'+star+'_'+self.model+self.elem+'_gc'+extra+'.pdf')
   
        if type(self.elem_ab)!=type(None):
            print(self.elem,'abundace: '+str(round(self.elem_ab[np.argmin(self.elem_EW_chi2)],2))+\
                  ' and microturbulence: '+str(round(self.elem_mic[np.argmin(self.elem_EW_chi2)],1)))
        ### FIRST APPROACH TO THE ABUNDANCE
      
        # plotting abundances figure: points for the different abundances that we got + line for mean value and shadow for the other. 
        # Dictionary of abundances with minimum chi2 (for all micro) for the used lines. Line-names as keys.
        self.ab_chi2min = {}; self.weight_chi2min = {}; self.ab_ion = {}
        self.lines_dic = {} # Dictionary with the lines (in a list) of each ion
        MeanMic = []
        
            ######################## FIGURE DEFINITION  #######################
        fig_sum, ax_sum = plt.subplots()
        fw.plot_style(ax_sum,size=12)
        ax_sum.set_xscale("log"); ax_sum.set_ylim(min(self.ab_FWgrid)+12,max(self.ab_FWgrid)+12)
        fig_sum.suptitle(self.star_name+', '+model+'\n' '$v\sin(i)$='+str(vsini)+\
        r' [km/s], v$_{\rm mac}$='+str(vmac)+' [km/s]')
        fig_sum.supxlabel(r'$EW\,[{\rm m}\AA]$')
        fig_sum.supylabel(r'$\log(\frac{\rm '+ion[0]+'}{\\rm H})+12[{\\rm dex}]$',fontsize=15)
        for ab_grid_tick in list(self.ab_FWgrid):
            ax_sum.axhline(ab_grid_tick+12,lw=1,c='grey',alpha=0.8,zorder=100,ls=':')
            ######################## FIGURE DEFINITION  #######################
            
        ############################### WORK PER ION ##############################
        for j,ion in enumerate(self.elem_list.keys()):
            ion_state = romanToInt(ion.replace(elem,''))
            self.lines_dic[ion] = []
            ############# Dictionaries for variable storage definition ############
            # microturb, abundance, chi2, EW_obs = {},{},{},{}
            microturb, abundance, chi2 = {},{},{}
            mic_Allmic,ab_Allmic, chi2_Allmic = {},{},{}
            micro_abundance, micro_EW, micro_chi2weight = {},{},{}
            ############# Dictionaries for variable storage definition ############
                  
            # If there are results from the line analysis in a given ion
            if len(self.elem_list[ion].keys())>0: 
                
                ######################## FIGURE DEFINITION  #######################
                
                        ############ General figure staff ############
                fig, axs = plt.subplots(len(self.elem_list[ion].keys())+1,sharex=True)
                fig.suptitle(self.star_name+', SpC='+SpC+'\n $v\sin(i)$='+str(vsini)+\
                r' [km/s], v$_{\rm mac}$='+str(vmac)+' [km/s] \n '+ion+', '+model)
                fig.subplots_adjust(right=0.8)
                fig.supxlabel(r'$\log(\frac{\rm '+ion[0]+'}{\\rm H})+12\,[{\\rm dex}]$',fontsize=15)
                fig.supylabel(r'$\xi_{\rm t}$ [km/s]', fontsize=18)
                [ax.set_xlim(min(self.ab_FWgrid)+12,max(self.ab_FWgrid)+12) for ax in axs]; [ax.set_ylim(float(min(self.mic_FWgrid).replace('VT','')),float(max(self.mic_FWgrid).replace('VT',''))) for ax in axs]
                
                fig_mic_ab,ax_mic_ab = plt.subplots(); fig_mic_ab.set_size_inches(9,7)
                fig_mic_ab.supxlabel(r'$\log(\frac{\rm '+ion[0]+'}{\\rm H})+12\,[{\\rm dex}]$',fontsize=15)
                fig_mic_ab.supylabel(r'$\xi_{\rm t}$ [km/s]',fontsize=18)
                ax_mic_ab.set_xlim(min(self.ab_FWgrid)+12,max(self.ab_FWgrid)+12);ax_mic_ab.set_ylim(float(min(self.mic_FWgrid).replace('VT','')),float(max(self.mic_FWgrid).replace('VT',''))) 
                fig_mic_ab.subplots_adjust(top=0.865)
                
                fw.plot_style([ax_mic_ab,*axs],size=12)          
                        ############ General figure staff ############
                
                # Grid with the FW grid steps marked.
                for ax in [ax_mic_ab,*axs]:
                    for mic_grid_tick in list(self.mic_FWgrid):
                        ax.axhline(float(mic_grid_tick.replace('VT','')),lw=1,c='grey',alpha=0.8,zorder=100,ls=':')
                    for ab_grid_tick in list(self.ab_FWgrid):
                        ax.axvline(ab_grid_tick+12,lw=1,c='grey',alpha=0.8,zorder=100,ls=':')
        
                ######################## FIGURE DEFINITION  #######################
                
                microturb[ion],abundance[ion],chi2[ion] = 0,0,0 # Starting mic, abun, and chi2 values for the total ion. 
                
                ####################### WORK PER LINE IN ION ######################
                for i,l in enumerate(self.elem_list[ion].keys()): 
                    line_EW = float(observed_EW_dic['EW_obs'][observed_EW_dic['line_name'] == l])
                    self.lines_dic[ion].append(l) # Including line studied for the given ion
                    #EW_obs[l] = observed_EW_dic['EW_obs'][l] # Including the observed EW of the studied line in a deicated dictionary 
                    ############### WORKING IN LINES OF A GIVEN ION ###############
                    line_mic = self.elem_list[ion][l]['microturbulences']
                    line_ab = self.elem_list[ion][l]['abundances']
                    line_chi2 = self.elem_list[ion][l]['chi2_EW']
                        ##### Normalizing chi2 to 1 and saving best 10% #####
                    max_level = sorted(line_chi2)[len(line_chi2)//10+1]
                    last_val = sorted(line_chi2)[len(line_chi2)//10]/max_level   
                    line_chi2 /= max_level
                    line_chi2[line_chi2>1] = 1
                        ##### Normalizing chi2 to 1 and saving best 10% ##### 
                        
                            #### Reshaping from linear to squared ####
                    x_size,y_size = len(np.unique(line_mic)),len(np.unique(line_ab))
                    mic_square = line_mic.reshape(x_size,y_size)
                    ab_square = line_ab.reshape(x_size,y_size)
                    chi2_square = line_chi2.reshape(x_size,y_size)    
                            #### Reshaping from linear to squared ####      
                    
                    ### Filling the dictionary (for line and ion)
                    microturb[l]= mic_square; abundance[l] = ab_square; chi2[l] = chi2_square
                    microturb[ion] += mic_square
                    abundance[ion] += ab_square
                    chi2[ion]      += chi2_square * self.line_weight[l]
                    # chi2 is already wighted by EW_obs + relative error. We weighted by line weight (to avoid overcounting transitions).REMOVED *observed_EW_dic['EW_obs'][l]                
                    
                    ##### Saving best parameters for ALL micros #####
                    mic_Allmic[l],ab_Allmic[l],chi2_Allmic[l] = [],[],[]
                    for k,min_mic in enumerate(mic_square):
                        min_mic = float(np.unique(min_mic))
                        min_chi2 = min(chi2_square[k]) # Finding the minimum chi2 for this mic
                        min_ab = ab_square[k][chi2_square[k]==min_chi2] # Corresponding abundance for the min(chi2) at a given micro.      
                        try: min_ab = float(min_ab)
                        except: min_ab = float(np.mean(min_ab))
                        ##### Best combinations for all micros #####
                        mic_Allmic[l].append(min_mic)
                        ab_Allmic[l].append(min_ab)  
                        chi2_Allmic[l].append(min_chi2)     
                        
                        ##### Saving ab and EW as a functions of micro #####
                        if str(min_mic) not in micro_abundance.keys():
                            micro_abundance[str(min_mic)],micro_EW[str(min_mic)],micro_chi2weight[str(min_mic)] = [],[],[]
                        micro_abundance[str(min_mic)].append(min_ab)   
                        micro_EW[str(min_mic)].append(line_EW) 
                        micro_chi2weight[str(min_mic)].append((1-min_chi2)*self.line_weight[l]) # Weight for averaging
                        #micro_chi2[str(min_mic)].append(1/min_chi2*line_weight[l]) # Weight for averaging
                    mic_Allmic[l], ab_Allmic[l], chi2_Allmic[l] = np.array(mic_Allmic[l]), np.array(ab_Allmic[l]), np.array(chi2_Allmic[l])
                    
                    #% 
                    ##### Plotting in the different figures #####
                    chi2_min_pos = np.argmin(chi2_Allmic[l]) # Mask of minimum chi2
                    
                    last_val = np.max(chi2[l][chi2[l]<1])
                    levels = np.linspace(np.min(chi2[l]),last_val, 21)
                    contourf_el = axs[i].contourf(abundance[l],microturb[l],chi2[l],levels=levels,cmap = cm.YlOrBr_r)
                    axs[i].plot(ab_Allmic[l],mic_Allmic[l],'b') # Line of minimum ab-mic pair for each mic.
                    axs[i].plot(ab_Allmic[l][chi2_min_pos],mic_Allmic[l][chi2_min_pos],'xr') # ab-mic of minimum chi2.
                    
                    ax_mic_ab.plot(ab_Allmic[l],mic_Allmic[l],lw=1,label=l+'\n EW='+str(line_EW)+r"m$\AA$") # Line of minimum ab-mic pair for each mic.
                    ax_mic_ab.plot(ab_Allmic[l][chi2_min_pos],mic_Allmic[l][chi2_min_pos],\
                            'X',markersize=6,markeredgewidth=0.5,c=colors[i])   
        
                    
                    axs[i].set_title(l+'\n EW='+str(line_EW)+r"m$\AA$",x=0.9,y=0.05)
                    ##### Plotting in the different figures #####
                ####################### WORK PER LINE IN ION ######################
                ###################### Saving ion parameters ######################
                microturb[ion] /= len(self.lines_dic[ion])   
                abundance[ion] /= len(self.lines_dic[ion])  
                chi2[ion] /= np.max(chi2[ion]) # Normalized to 1
                ###################### Saving ion parameters ######################
        
                ####### Saving best combinations for all micros for the ion #######
                mic_Allmic[ion],ab_Allmic[ion],chi2_Allmic[ion] = [],[],[]
                for k,min_mic in enumerate(microturb[ion]):
                    min_mic = float(np.unique(min_mic))
                    min_chi2 = min(chi2[ion][k]) # Finding the minimum chi2 for this mic
                    min_ab = float(abundance[ion][chi2[ion]==min_chi2]) # Corresponding abundance for the min(chi2) at a given micro.      
                    
                    ##### Best combinations for all micros #####
                    mic_Allmic[ion].append(min_mic)
                    ab_Allmic[ion].append(min_ab)  
                    chi2_Allmic[ion].append(min_chi2)   
                mic_Allmic[ion], ab_Allmic[ion], chi2_Allmic[ion] = np.array(mic_Allmic[ion]), np.array(ab_Allmic[ion]), np.array(chi2_Allmic[ion])  
                micro_abundance = np.transpose(np.array(list(micro_abundance.values())))
                micro_chi2weight = np.transpose(np.array(list(micro_chi2weight.values())))
                micro_EW = np.transpose(np.array(list(micro_EW.values())))
                ####### Saving best combinations for all micros for the ion #######
        
                            ### Plotting total map for a given ion ###
                i += 1 
                chi2_min_pos = np.argmin(chi2_Allmic[ion])
                last_val = np.max(chi2[ion][chi2[ion]<1])
                levels = np.linspace(np.min(chi2[ion]),last_val, 21)
                contourf_el = axs[i].contourf(abundance[ion],microturb[ion],chi2[ion],levels=levels,cmap = cm.YlOrBr_r) 
                cbar_ax = fig.add_axes([0.85, 0.15, 0.03, 0.7])
                cbar = fig.colorbar(contourf_el, cax=cbar_ax, ticks=[0,last_val])
                cbar.ax.set_yticklabels(['0', 'best 10%']) 
                cbar.set_label(r'$\chi^2=\left(\frac{EW_{\rm mod}-EW_{\rm obs}}{EW_{\rm obs}}\right)^2\Delta EW$',\
                               loc='center', fontsize=14,labelpad=-20)
                ############ General figure staff ############
                
                axs[i].plot(ab_Allmic[ion],mic_Allmic[ion],'b') # Line of minimum ab-mic pair for each mic.
                axs[i].plot(ab_Allmic[ion][chi2_min_pos],mic_Allmic[ion][chi2_min_pos],'xr') # ab-mic of minimum chi2.
                            ### Plotting total map for a given ion ###
                # Measuring the weighted standard deviation of all lines 
                mean_ab = np.average(micro_abundance,axis=0,weights=micro_chi2weight)
                var_ab = (micro_abundance-mean_ab)**2
                variance_ab = np.average(var_ab,axis=0,weights=micro_chi2weight)
                var_ab = np.sqrt(variance_ab)
                
                #x1, x2 = ab_Allmic[ion] - var_ab, ab_Allmic[ion] + var_ab # Confidence limits from map minimum 
                x1, x2 = mean_ab - var_ab, mean_ab + var_ab # Confidence limits from weighted average
                
                i_min_sct = np.argmin(var_ab) # Minimazing microturbulence scatter
                
                ax_mic_ab.plot(ab_Allmic[ion],mic_Allmic[ion],'--k',lw=2.5,\
                               label=r'Map min. pairs') # Minimum from total map
                ax_mic_ab.plot(mean_ab,mic_Allmic[ion],':k',lw=2.5,label='Weighted min. pairs') # Minimum from weighted average
                ax_mic_ab.fill_betweenx(mic_Allmic[ion],x1,x2,color='k',alpha=0.2,label='Confidence area') # Confidence area (from weighted average)
                #ax_mic_ab.fill_betweenx(mic_Allmic[ion],x1,x2,color='k',alpha=0.2) # Confidence area (from map minimum)
        
                if type(ab_Allmic[ion][i_min_sct])==np.float64:
                    ax_mic_ab.plot(ab_Allmic[ion][i_min_sct],mic_Allmic[ion][i_min_sct],\
                            '^k',markersize=6,markeredgewidth=1.5) # Minimum from minimizing scattering
                    ax_mic_ab.plot(mean_ab[i_min_sct],mic_Allmic[ion][i_min_sct],\
                            '^k',markersize=8,markeredgewidth=1.5) # Minimum from minimizing scattering
                else: # If there is more than one minimum, we choose the first one. 
                    ax_mic_ab.plot(ab_Allmic[ion][i_min_sct][0],mic_Allmic[ion][i_min_sct][0],\
                            '^k',markersize=6,markeredgewidth=1.5) # Minimum from minimizing scattering
                    ax_mic_ab.plot(mean_ab[i_min_sct][0],mic_Allmic[ion][i_min_sct][0],\
                            '^k',markersize=8,markeredgewidth=1.5) # Minimum from minimizing scattering
                
                i_min_chi2 = np.argmin(chi2_Allmic[ion]) # Minimizing chi2 in the total map.
        
                ax_mic_ab.plot(ab_Allmic[ion][i_min_chi2],mic_Allmic[ion][i_min_chi2],\
                        'xk',markersize=8,markeredgewidth=1.5) # Minimum from minimizing chi2
                
                # Cleaning: Refuse lines that are beyond 3 times the scatter in the minimum point (sigma_clipping)
                old_lenght = np.shape(micro_abundance)[0]
                new_lenght = old_lenght + 1
                var_ab_sc = np.copy(var_ab)
                while new_lenght > old_lenght:
                    i_min_sct = np.argmin(var_ab_sc) # Minimazing microturbulence scatter
                    ab_min = ab_Allmic[ion][i_min_sct] - 3 * np.nanmin(var_ab_sc) 
                    ab_max = ab_Allmic[ion][i_min_sct] + 3 *  np.nanmin(var_ab_sc) 
                    # Creating the mask with the new limits
                    mask_sc = ((micro_abundance[:,i_min_sct] >= ab_min) \
                         & (micro_abundance[:,i_min_sct] <= ab_max))
                    # Passing the mask to an array with the corresponding shape
                    mask_sc_array = np.transpose(np.array([mask_sc for k in range(np.shape(micro_abundance)[1])]))
                                    
                    # New mean + sct values. We use np.where to conserve the dimession.
                    micro_abundance_sc = np.where(mask_sc_array, micro_abundance, np.nan)
                    micro_chi2weight_sc = np.where(mask_sc_array, micro_chi2weight, np.nan)
                    mean_ab_sc = np.average(micro_abundance_sc,axis=0,weights=micro_chi2weight_sc)
                    var_ab_sc = (micro_abundance_sc - mean_ab_sc)**2
                    variance_ab_sc = np.average(var_ab_sc,axis=0,weights=micro_chi2weight_sc)
                    var_ab_sc = np.sqrt(variance_ab_sc)
                    # Givning new values to the lengths
                    old_lenght = new_lenght; new_lenght = mask_sc.sum()
                i_min_sct = np.argmin(var_ab_sc)
                x1_sc, x2_sc = mean_ab_sc - var_ab_sc, mean_ab_sc + var_ab_sc 
                if max(abs(x1_sc-x1))+max(abs(x2_sc-x2)) >= 1e-3:
                    ax_mic_ab.plot(mean_ab_sc,mic_Allmic[ion],'k',lw=2.5,label='Clipped weighted minimum')
                    ax_mic_ab.fill_betweenx(mic_Allmic[ion],x1_sc,x2_sc,color='b',alpha=0.2,label='Clipped confidence area')
                    if len(mean_ab[i_min_sct])==1:
                        ax_mic_ab.plot(mean_ab[i_min_sct],mic_Allmic[ion][i_min_sct],'^k',\
                            markersize=8,markeredgewidth=1.5)
                    else:
                        ax_mic_ab.plot(mean_ab[i_min_sct][0],mic_Allmic[ion][i_min_sct][0],'^k',\
                            markersize=8,markeredgewidth=1.5)
                
                # esquina superior derecha en nombre de la estrella, el SpC, la vsini, vmac, y el nombre del modelo 
                # Writting parameters in the upper right part 
                fig_mic_ab.suptitle(self.star_name+', SpC='+SpC+'\n $v\sin(i)$='+str(vsini)+\
                r' [km/s], v$_{\rm mac}$='+str(vmac)+' [km/s] \n '+\
                                   ion+', '+model,fontsize=12)
                
                            ### Legend definition & save fig ###
                
                ax_mic_ab.plot([],[],'^k',markersize=8,markeredgewidth=1.5,label=r'Minimum $\sigma$'+' \n '+\
                    '$\epsilon_{\\rm '+ion+'}='+str(round(mean_ab_sc[i_min_sct],2))+'\pm'+\
                     str(round(var_ab_sc[i_min_sct],2))+r'$; $\xi ='+str(mic_Allmic[ion][i_min_sct])+'$')
                    
                ax_mic_ab.plot([],[],'X',color='grey',markersize=7,markeredgewidth=0.5,label=r'Line minimum $\chi^2$'+' \n '+\
                    '$\epsilon_{\\rm '+ion+'}='+str(round(ab_Allmic[ion][i_min_chi2],2))+'\pm'+\
                     str(round(var_ab_sc[i_min_chi2],2))+r'$; $\xi ='+str(mic_Allmic[ion][i_min_chi2])+'$')
                
                ax_mic_ab.legend()
                if save_im:
                    fig_mic_ab.savefig(save_im_path+star+'_'+model+'_'+ion+'_micro_ab'+extra+'.pdf')
                    fig.savefig(save_im_path+star+'_'+ion+'_'+model+extra+'.pdf')
                            ### Legend definition & save fig ###
                
                             ### Fixing microturbulence for the result ###
                # Once we know the commun microturbulence we like, we get the correspondent abundance. 
                if scatter_minimum:
                    i_result = i_min_sct
                else:
                    i_result = i_min_chi2
                micro_min_sct = mic_Allmic[ion][i_result] # Microturbulences mic
                star_analisys.micro_func()
                micL = round((star_analisys.mic_L))
                micsct = (micro_min_sct)
                i_micL = mic_Allmic[ion] == micL; i_micsct = mic_Allmic[ion] == micsct
                ab_ion_L, std_ab_ion_L = round(mean_ab_sc[i_micL][0],2), round(var_ab_sc[i_micL][0],2)
                ab_ion_sct, std_ab_ion_sct = round(mean_ab_sc[i_micsct][0],2), round(var_ab_sc[i_micsct][0],2)
                
                if abs(ab_ion_L-ab_ion_sct) >= np.sqrt(std_ab_ion_L**2 + std_ab_ion_sct**2):
                    print('\033[41mWARNING! \033[0m'+'\033[33mDifference between abundances with L_micro vs \
                          scatter_micro greater than uncertainty \033[0m')
                self.mic_min_sct_ion[ion] = micro_min_sct
                # Trying to fix micro with multiplet minimum scatter  
                for common_transition in list(tran_dic.keys()):
                    lines_common_transition = tran_dic[common_transition]
                    if len(lines_common_transition) >= 3:
                        multi_ab,multi_chi2 = [],[]
                        multi_study = False
                        for lin in lines_common_transition:
                            if lin in list(mic_Allmic.keys()):
                                multi_study = True
                                multi_ab.append(ab_Allmic[lin])
                                multi_chi2.append((1-np.array(chi2_Allmic[lin]))*self.line_weight[lin])
                        if multi_study:
                            multi_ab, multi_chi2 = np.array(multi_ab), np.array(multi_chi2)
                            multi_mean_ab = np.average(multi_ab,axis=0,weights=multi_chi2)
                            multi_var_ab = (multi_ab-multi_mean_ab)**2
                            multi_variance_ab = np.average(multi_var_ab,axis=0,weights=multi_chi2)
                            multi_var_ab = np.sqrt(multi_variance_ab)
                            i_min_sct_multi = np.argmin(multi_var_ab)
                            self.mic_min_sct_multiplet[common_transition] = mic_Allmic[ion][i_min_sct_multi]
                
                if fix_micro:
                    if type(fix_micro) == dict:
                        mic = fix_micro[star]
                    else:
                        mic = fix_micro
                        
                elif mic_L:
                    print('Fixing micro according to spectral luminosity')
                    mic = micL
                    
                else:
                    try: mic = round(np.mean(list(self.mic_min_sct_multiplet.values())))
                    except: mic = micsct
                    
                print('Micro is',mic)
                MeanMic.append(mic) # Saving microturbulence value for the studied ion.
                
                            #### Saving results in a  dictionary. ####

                # Save all abundance informatio in a dictionary (we can recompute different abundance approaches with this dictionary)
                for lin,ab,use,w in zip(np.array(list(self.elem_list[ion].keys())),micro_abundance,mask_sc,micro_chi2weight):    
                    self.ab_chi2min[lin] = ab
                    if use: 
                        self.weight_chi2min[lin] = w
                    else: self.weight_chi2min[lin] = -1*w

                self.ab_ion[ion] = {'ab': mean_ab_sc, 'uncertainty': var_ab_sc}
        self.micro_array =  mic_Allmic[ion]       
        ############################### WORK PER ION ##############################
        
                 ####### Ploting results in the summary figure #######
        
        # N: Number of lines for a given ion
        # ion_mean: Weighted average abundance for the given ion
        # ion_std: Error in abundance value for the given ion
        i_result = mic_Allmic[ion] == mic # Recomputing the mask with the corresponding microturbulence
        N = mask_sc.sum()

        for ion in list(self.ab_ion.keys()):
            ion_state = romanToInt(ion.replace(elem,''))
            ax_sum.axhspan([self.ab_ion[ion]['ab']-self.ab_ion[ion]['uncertainty']][0][i_result][0],\
                           [self.ab_ion[ion]['ab']+self.ab_ion[ion]['uncertainty']][0][i_result][0], facecolor=colors[ion_state], alpha=0.2)
            
            sum_ab, sum_n = 0,0
            for lin in list(self.elem_list[ion].keys()):
                line_EW = float(observed_EW_dic['EW_obs'][observed_EW_dic['line_name'] == lin])
                ax_sum.annotate(lin,xy=(line_EW,self.ab_chi2min[lin][i_result][0]),fontsize=8,color=colors[ion_state])  
                if sum(self.weight_chi2min[lin]) > 0:
                    # Lines used for the abundance determination 
                    ax_sum.scatter(line_EW,self.ab_chi2min[lin][i_result][0],c=colors[ion_state],\
                        s=mic)
                    sum_ab += self.ab_chi2min[lin][i_result][0]
                    sum_n += 1
                else:
                    #Lines removed during sigma-clipping for the abundance determination  
                    ax_sum.scatter(line_EW,self.ab_chi2min[lin][i_result][0],\
                        micro_abundance[:,i_result][np.invert(mask_sc)],\
                        facecolors='none', edgecolors=colors[ion_state],\
                        s=mic)  
                        #,label=ion+', N='+str(n)+r', $\epsilon=$'+str(ion_mean_ab)+r'$\pm$'+str(ion_std_ab)+r', $\xi=$'+mic+'km/s')
            
            ax_sum.scatter([],[],c=colors[ion_state],s=mic,label=ion+', N='+str(N)+r', $\epsilon=$'+\
                str(round(self.ab_ion[ion]['ab'][i_result][0],2))+r'$\pm$'+\
                str(round(self.ab_ion[ion]['uncertainty'][i_result][0],2))+\
                ' \n'+r'$\xi=$'+str(mic)+'km/s')
            ax_sum.axhline(self.ab_ion[ion]['ab'][i_result][0],lw = 2, color=colors[ion_state],label = 'Weighted average')
            ax_sum.axhline(sum_ab/sum_n,ls=':',color=colors[ion_state],label = 'Mean values \n WITHOUT weighted')
                     ####### Ploting results in the summary figure #######        
        
        # Plotting lines with EW NOT INCLUDED IN THE STUDY (from que scratch)
        for l in range(len(observed_EW_dic['line_name'])):
            if observed_EW_dic['Activation'][l]==0 and observed_EW_dic['EW_obs'][l]>=10: 
                print('IMPORTANT: ',mic)
                lin = observed_EW_dic['line_name'][l]
                EW_obs = observed_EW_dic['EW_obs'][l]#observed_EW_dic['EW_obs'][l]
                ion = lin[:re.search(r"\d",lin).start()]
                ion_state = romanToInt(ion.replace(elem,''))
                model_line = {'abundance' : [], 'microturbulence' : [], 'EW' : []}
                ################ Generating the original model dictionary ################
                for ab in self.FW_model.lines[lin].abundances.keys(): # Computing pseudo-xi
                    for mic_for in self.FW_model.lines[lin].abundances[ab].microturbulences.keys():
                        try: 
                            model_line['EW'].append(self.FW_model.get_profile(lin,ab,mic_for).EW)#*1e6#/EW_obs)) # We multiply by 1e6 due to the error when saving EW
                            model_line['abundance'].append(float(ab)+12)
                            model_line['microturbulence'].append(float(mic_for[-3:]))
                        except: pass
                ab_FW,mic_FW,ew_FW = np.array(model_line['abundance']),np.array(model_line['microturbulence']),np.array(model_line['EW'])
                ########## Interpolation ##########
                ab_interp = np.linspace(np.nanmin(ab_FW),np.nanmax(ab_FW),step_number_abundance)
                mic_interp = np.linspace(np.nanmin(mic_FW),np.nanmax(mic_FW),step_number_microturbulence)
             
                ab_square,mic_square = np.meshgrid(ab_interp,mic_interp)
                EW = griddata((ab_FW,mic_FW), ew_FW, (ab_square,mic_square), method='cubic')
                chi2_square = ((EW_obs-EW)/EW_obs)**2
                ab_linearized = ab_square.reshape(len(ab_interp)*len(mic_interp))
                mic_linearized = mic_square.reshape(len(ab_interp)*len(mic_interp))
                chi2_linearized = chi2_square.reshape(len(ab_interp)*len(mic_interp))
                mic_ter = sgf.find_nearest(mic_linearized,mic)
                resulting_abund = ab_linearized[mic_linearized==mic_ter][np.argmin(chi2_linearized[mic_linearized==mic_ter])]
                ########## Interpolation ##########
                
                ax_sum.plot(EW_obs,resulting_abund,'x',c=colors[ion_state])
                ax_sum.annotate(lin,xy=(EW_obs,resulting_abund),fontsize=8,color=colors[ion_state])
        
        ab_chi2min, weight_chi2min = [], []
        for ab, weight in zip(list(self.ab_chi2min.values()),list(self.weight_chi2min.values())):
            if sum(weight) > 0:
                ab_chi2min.append(ab), weight_chi2min.append(weight)
        ab_chi2min, weight_chi2min = np.array(ab_chi2min), np.array(weight_chi2min)
        self.weighted_abundance = np.around(np.average(ab_chi2min,weights=weight_chi2min,axis=0),2)
        # Average abundance using ALL the lines. 
        variance = np.average((ab_chi2min-self.weighted_abundance)**2,\
          weights=weight_chi2min,axis = 0)
        self.weighted_uncertainty = np.around(np.sqrt(variance),2)#round(np.std(ab_tot),2)
        
                        ### Legend definition & save fig ###
        ax_sum.plot([],[],' ',label=r'$\epsilon='+str(self.weighted_abundance[self.micro_array == mic][0])+\
        '\pm'+str(self.weighted_uncertainty[self.micro_array == mic][0])+'$')
        #ax_sum.plot([],[],' ',label=r'$\xi_t='+str(np.nanmean(MeanMic))+'$ km/s')        
        box = ax_sum.get_position()
        ax_sum.set_position([box.x0, box.y0, box.width * 0.75, box.height])
        ax_sum.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        if save_im: fig_sum.savefig(save_im_path+star+'_'+model+'_EW_abundance_gc_fabrice.pdf')
                        ### Legend definition & save fig ###
                        
        
        return {'Micro': mic,\
                'Weighted abundance': self.weighted_abundance[self.micro_array == mic][0],\
                'Weighted uncertainty': self.weighted_uncertainty[self.micro_array == mic][0]}

    
if __name__ is '__main__':

    '''
    star = 'HD55879'
    model = 'T310g350He10Q1350b10'
    SpC = 'O9.7III'; Teff = 31.5; logg = 3.54
    vsini = 28;vmac = 61
    
    '''
    star = 'HD207538'
    model = 'T330g380He10Q1350b10'
    SpC = 'O9.7IV'; Teff = 31.9; logg = 3.81
    vsini = 29;vmac = 50
    
    elem = 'N'
    
    
    
    ### Measured EW in observations ###
    read_EW = pd.read_csv('/home/cms/Desktop/PhD/0_PROJECT/MONTPELLIER/STAR_ANALYSIS/'+star+'/'+star+'_N_gc_fabrice.csv')
    ### Measured EW in observations ###
    #%
    
    
    star_analisys = abundance_determination(star,Teff,logg)
    results = star_analisys.abundance_determination_func('N',model,vsini,vmac,SpC=SpC,observed_EW_dic=read_EW)


    
