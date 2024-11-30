# -*- coding: utf-8 -*-
"""
@author: C. Martinez Sebastian
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import formulas as fmr
#import tools as tls
import pickle
from matplotlib.lines import Line2D
import spec as sp
import line_broadening as bd
from db import *

def plot_style(AX,size=18):
    if type(AX) == list:
        for ax in AX: 
            ax.tick_params(which='minor', length=2, direction='in', top='on', right='on')
            ax.tick_params(which='major', length=5, direction='in', top='on', right='on')
            ax.xaxis.set_minor_locator(AutoMinorLocator(5))
            ax.yaxis.set_minor_locator(AutoMinorLocator(5))
            ax.xaxis.set_tick_params(labelsize=size)
            ax.yaxis.set_tick_params(labelsize=size)
    else: 
        ax = AX      
        ax.tick_params(which='minor', length=2, direction='in', top='on', right='on')
        ax.tick_params(which='major', length=5, direction='in', top='on', right='on')
        ax.xaxis.set_minor_locator(AutoMinorLocator(5))
        ax.yaxis.set_minor_locator(AutoMinorLocator(5))
        ax.xaxis.set_tick_params(labelsize=size)
        ax.yaxis.set_tick_params(labelsize=size)
        
def considered_transitions():
    '''
    Parameters
    ----------
    NONE

    Returns
    -------
    transitions : TYPE
        DESCRIPTION.
    central_wavelenths : TYPE
        DESCRIPTION.

    '''
    ############### Considered transition ############### 
    with open(dir().modelgrid+'non_repeated_lines_new.txt') as f:
        lines = f.readlines()
    f.close()
    considered_lines = []
    transitions = {}
    for l in lines:
        try: 
            if int(l.split()[1].replace(",",""))==1:
                considered_lines.append(l.split()[0].replace(",",""))
        except: pass
        if "#" not in l:
            try:
                transition_levels = []
                transition_levels.append(l.split()[1].replace(",",""))
                for i in range(int((len(l.split())-4)/4)):
                    transition_levels.append(l.split()[4*i+4].replace(",","")+'-'+\
                                             l.split()[4*i+5].replace(",","")+'-'+\
                                             l.split()[4*i+6].replace(",",""))
                transitions[l.split()[0].replace(",","")] = transition_levels
            except: pass
    ############### Considered transition ############### 
    
    ############### Central_wavelenth for each transition ############### 
    DATA = pd.read_csv(dir().modelgrid+'LINES_CNOSi_new_coll.dat',sep='\s+',names=\
                               ["i_l","f_l","fine","CW","oscillator_strengths",\
                                "radiative_damping","c7","c8","collisional_damping",\
                                "sum_gi"],on_bad_lines='skip')
    DATA_H = pd.read_csv(dir().modelgrid+'LINES_HHe_new_coll.dat',sep='\s+',names=\
                               ["i_l","f_l","fine","CW"],on_bad_lines='skip')#,error_bad_lines=False)
    central_wavelenths = {}
    for i, il in enumerate(DATA.i_l):
        if ':T' in il:
            pass
        else:
            try: central_wavelenths[il+'-'+DATA.f_l[i]+'-'+DATA.fine[i]] = DATA.CW[i]
            except: pass
    for i, il in enumerate(DATA_H.i_l):
        if ':T' in il:
            pass
        else:
            #print(il)
            try: central_wavelenths[il+'-'+DATA_H.f_l[i]+'-'+DATA_H.fine[i]] = DATA_H.CW[i]
            except: pass
    ############### Central_wavelenth for each transition ############### 
  
    return transitions, central_wavelenths

def plot_spectra(XY,transitions, central_wavelenths, wv1 = 3800,wv2 = 7000):
    '''
    

    Parameters
    ----------
    XY : list,
        All the components must have the same dimesion
    transitions : TYPE
        DESCRIPTION.
    central_wavelenths : TYPE
        DESCRIPTION.
    wv1 : TYPE, optional
        DESCRIPTION. The default is 3800.
    wv2 : TYPE, optional
        DESCRIPTION. The default is 7000.

    Returns
    -------
    None.

    '''
    #plt.close('all')
    fig, ax = plt.subplots(len(XY)+1,sharex=True,sharey=True)
    fig.supxlabel(r'Wavelength [$\AA$]'),fig.supylabel(r'Normalized flux [ad]')
    ax[0].xaxis.set_minor_locator(AutoMinorLocator(5))
    ax[0].yaxis.set_minor_locator(AutoMinorLocator(4))
    xy_tot = np.zeros(len(XY[0][1]))
    for i,xy in enumerate(XY):
        ax[i].tick_params(which='both', direction='in',bottom=True,top=True,left=True,right=True)
        ax[i].plot(xy[0],xy[1]+1),ax[i].set_ylabel(xy[2])
        globals()[xy[2][0]] = ax[i]
        xy_tot += xy[1]
    ax[i+1].tick_params(which='both', direction='in',bottom=True,top=True,left=True,right=True)
    ax[i+1].plot(xy[0],xy_tot+1),ax[i].set_ylabel('Total')      
    
    ############### Where should we print the transitions ############### 
    elem = []
    for k in transitions.keys():
        if k[0] not in elem:
            elem.append(k[0])
            globals()[k[0]+'_not_considered'],globals()[k[0]+'_considered'] = [],[]
        if transitions[k][0]=='1':
            globals()[k[0]+'_considered'].extend(transitions[k][1:])
        else:
            globals()[k[0]+'_not_considered'].extend(transitions[k][1:])
    for E in elem:
        globals()[E+'_not_considered_def'] = []
        globals()[E+'_not_considered'] = list(dict.fromkeys(globals()[E+'_not_considered']))
        globals()[E+'_considered'] = list(dict.fromkeys(globals()[E+'_considered']))
        for t in globals()[E+'_not_considered']:
            if t not in globals()[E+'_considered']:
                globals()[E+'_not_considered_def'].append(t)
    ############### Where should we print the transitions ############### 
    for E in elem:
        for t in globals()[E+'_considered']:
            try: 
                globals()[E].text(float(central_wavelenths[t])-0.5,1.1,t,c='b',rotation = 70, fontsize = 'x-small', clip_on=True)
                globals()[E].vlines(float(central_wavelenths[t]),1.05,1.09,colors='b')
            except: print(t)
        for t in globals()[E+'_not_considered_def']:
            try: 
                globals()[E].text(float(central_wavelenths[t])-0.5,1.1,t,c='r',rotation = 70, fontsize = 'x-small', clip_on=True)
                globals()[E].vlines(float(central_wavelenths[t]),1.05,1.09,colors='r')
            except: print(t)
    ax[0].set_xlim(wv1,wv2)
    return xy_tot

def abundance(E,E_sun):
    A = E_sun-0.01,E_sun,E_sun+0.01
    N = fmr.find_nearest(A,E)
    if N[1]==0: z = 'm010'
    elif N[1]==1: z = 'p000'
    elif N[1]==2: z = 'p010'
    return z

def model_comparation(star_spectrum,E_C,E_N,E_O,model,vmic,vsini,vmac,R,beta=1.5,rv0=0,plot=False):
    pkl_folder=dir().modelgrid+dir().modelgrid.split('/')[-2]+'_pkl'
    star = sp.spec(star_spectrum,rv0=rv0) # Radial velocity from SIMBAD.
    obs_wave, obs_flux = star.wave, star.flux
    z_C,z_N,z_O = abundance(E_C,-3.70+12),abundance(E_N,-4.40+12),abundance(E_O,-3.25+12)
    # We have to consider the same abundance of p... in N and in He+H. 
    with open(pkl_folder+model[:-3]+'_C'+z_C+'.pkl', 'rb') as inp:   
        C_lines = pickle.load(inp)
        inp.close()
    with open(pkl_folder+model[:-3]+'_N'+z_N+'.pkl', 'rb') as inp:   
        N_lines = pickle.load(inp)
        inp.close()
    with open(pkl_folder+model[:-3]+'_O'+z_O+'.pkl', 'rb') as inp:   
        O_lines = pickle.load(inp)
        inp.close()
    with open(pkl_folder+model[:-3]+'_H_He_'+z_N+'.pkl', 'rb') as inp:   
        H_lines = pickle.load(inp)
        inp.close()
    if type(vmic)!=str:
        microturb = H_lines.microturb.__dict__.keys()
        m, vt = [],[]
        for mic in microturb: 
            m.append(float(mic[2:])), vt.append(mic)
        vmic = vt[fmr.find_nearest(m,vmic)[1]]
        print(vmic)
    # To see the possible microturbulences: C_lines.microturb.__dict__.keys()
    
    C_L = getattr(C_lines.microturb,vmic)
    N_L = getattr(N_lines.microturb,vmic)
    O_L = getattr(O_lines.microturb,vmic)
    H_L = getattr(H_lines.microturb,vmic)
    #C_L = C_lines.microturb.VT007
    #N_L = N_lines.microturb.VT007
    #O_L = O_lines.microturb.VT007
    #H_L = H_lines.microturb.VT007
    
    transitions, central_wavelenths = considered_transitions()    
    
    wv_C,C_flux,E_C = C_L['total']['wavelength'],C_L['total']['flux'],C_lines.header.E_C
    wv_N,N_flux,E_N = N_L['total']['wavelength'],N_L['total']['flux'],N_lines.header.E_N
    wv_O,O_flux,E_O = O_L['total']['wavelength'],O_L['total']['flux'],O_lines.header.E_O
    wv_H,H_He_flux,Y_He = H_L['total']['wavelength'],H_L['total']['flux'],H_lines.header.Y_He
    #T_flux = C_flux + N_flux + O_flux + H_He_flux
    
    XY = [wv_C,C_flux,r'C, $\epsilon_{C}=$'+str(E_C)],\
        [wv_N,N_flux,r'N, $\epsilon_{N}=$'+str(E_N)],\
        [wv_O,O_flux,r'O, $\epsilon_N=$'+str(E_O)],\
        [wv_H,H_He_flux,r'H+He, $Y_{He}=$'+str(Y_He)]
    
    T_flux = plot_spectra(XY,transitions, central_wavelenths)
    # Convolution
    
    with open('line_broadening_divisions.txt','r') as f:
        lines = f.readlines()
    f.close()
    mask = []
    for l in lines:
        if '#' not in l:
            l.replace('\n','')
            mask.append([float(l.split()[0]),float(l.split()[1])])
    FL_conv = np.zeros(len(T_flux))
    for m in mask:
        l0 = np.mean(m)
        M = (wv_O>m[0]) & (wv_O<m[1])
        WAVE, FLUX = wv_O[M],T_flux[M]
        WV, FL = bd.final_resolution(WAVE, FLUX, l0, vsini, vmac, R, beta=beta)
        WV, FL = tls.interpol(WV, FL, np.mean(np.diff(WAVE)),start=round(WAVE[0],2),stop=round(WAVE[-1]+0.005,2))
        FL_conv[M] += FL

    wv,FL_conv = tls.interpol(wv_O,FL_conv,np.mean(np.diff(obs_wave)))
    XY_obs = obs_wave, obs_flux,{'color':'b'}
    XY_FW_conv = wv,FL_conv+1,{'color':'r'}
    XY = XY_obs,XY_FW_conv
    if plot:
        tls.plot(XY)
        plt.xlim(min(wv_O),6800)
    return XY

def fil_7(XY,title,color,transitions,vmax=1.1,vmin=0.6,models = None):
    XY_obs = XY[0]
    try: obs_wave,obs_flux,ot_1 = XY_obs
    except: obs_wave,obs_flux = XY_obs
    wv,FL_conv,ot = {},{},{}
    for i in range(1,len(XY)):
        try: wv[i],FL_conv[i],ot[i] = XY[i]
        except: wv[i],FL_conv[i] = XY[i]
    
    fig = plt.figure(figsize=(13,21))
    ax1,ax5,ax10,ax14,ax17,ax20,ax23 = plt.subplot(7,4,1),plt.subplot(7,5,6),plt.subplot(7,4,9),\
        plt.subplot(7,3,10),plt.subplot(7,3,13),plt.subplot(7,3,16),plt.subplot(7,2,13)
    ax2,ax3,ax4 = plt.subplot(7,4,2,sharey =ax1),plt.subplot(7,4,3,sharey =ax1),plt.subplot(7,4,4,sharey =ax1)
    ax6,ax7,ax8,ax9=  plt.subplot(7,5,7,sharey =ax5),plt.subplot(7,5,8,sharey =ax5),plt.subplot(7,5,9,sharey =ax5),plt.subplot(7,5,10,sharey =ax5)
    ax11,ax12,ax13 = plt.subplot(7,4,10,sharey =ax10),plt.subplot(7,4,11,sharey =ax10),plt.subplot(7,4,12,sharey =ax10)
    ax15,ax16 = plt.subplot(7,3,11,sharey =ax14),plt.subplot(7,3,12,sharey =ax14)
    ax18,ax19 = plt.subplot(7,3,14,sharey =ax17),plt.subplot(7,3,15,sharey =ax17)
    ax21,ax22 = plt.subplot(7,3,17,sharey =ax20),plt.subplot(7,3,18,sharey =ax20)
    ax24 = plt.subplot(7,2,14,sharey =ax23)
    
    
    xlim = [4079,4121],[4319,4361],[4839,4881],[6539,6581],\
        [4375,4405],[4910,4940],[4010,4040],[4455,4485],[4700,4730],\
        [4185,4215],[4525,4555],[4670,4700],[6670,6700],\
        [4498,4528],[6365,6395],[4596,4626],\
        [4252,4282],[5680,5710],[5786,5816],\
        [4400,4430],[3946,3976],[5577,5607],\
        [4620,4660],[4045,4090]
    lab = r'$H\delta+NIII4097$',r'$H\gamma$',r'$H\beta$',r'$H\alpha$',\
        r'$HeI4387$',r'$HeI4922$',r'$HeI+II4026$',r'$HeI4471$',r'$HeI4713$',\
        r'$HeI4200$',r'$HeII4541$',r'$HeII4686$',r'$HeI6678+HeII6683$',\
        r'$NIII4511/15$',r'$NIV6380$',r'$NV4603/19$',\
        r'$CII4267$',r'$CIII5695$',r'$CIV5801$',\
        r'$OII4414/16$',r'$OIII3961$',r'$OIII5592$',\
        r'$NII4630+NIII4634/40+CIII4647$',r'$NIV4058+OII4075$'
    for i in range(24):
        M_obs = (obs_wave>=xlim[i][0]) & (obs_wave<=xlim[i][1])
        if i==0 or i==4 or i==9 or i==13 or i==16 or i==19 or i==22:
            locals()['ax'+str(i+1)].tick_params(which='both', direction='in',bottom=True,top=True,left=True,right=True)
        else:
            locals()['ax'+str(i+1)].tick_params(which='both', direction='in',bottom=True,top=True,left=True,labelleft=False,right=True)
        locals()['ax'+str(i+1)].xaxis.set_minor_locator(AutoMinorLocator(5))
        locals()['ax'+str(i+1)].yaxis.set_minor_locator(AutoMinorLocator(4))
        locals()['ax'+str(i+1)].plot(obs_wave[M_obs],obs_flux[M_obs])
        for j in range(1,len(XY)):
            M_FW = (wv[j]>=xlim[i][0]) & (wv[j]<=xlim[i][1])
            locals()['ax'+str(i+1)].plot(np.array(wv[j])[M_FW],np.array(FL_conv[j])[M_FW])
            #print(M_FW)
        locals()['ax'+str(i+1)].set_title(lab[i])
        ad = np.diff(locals()['ax'+str(i+1)].set_ylim())
        locals()['ax'+str(i+1)].set_ylim(bottom=min(max(np.nanmin(obs_flux[M_obs]),vmin),np.nanmin(np.array(FL_conv[j])[M_FW]),locals()['ax'+str(i+1)].set_ylim()[0]),\
                                          top=max(1+0.25*ad,locals()['ax'+str(i+1)].set_ylim()[1],min(np.nanmax(np.array(FL_conv[j])[M_FW]),vmax)))
        if i==0 or i==4 or i==9 or i==13 or i==16 or i==19 or i==22:
            locals()['ax'+str(i+1)].set_ylim(bottom=min(max(np.nanmin(obs_flux[M_obs]),vmin),np.nanmin(np.array(FL_conv[j])[M_FW])),\
                                          top=max(min(np.nanmax(obs_flux[M_obs]),vmax,),np.nanmax(np.array(FL_conv[j])[M_FW]),1+0.2*ad))            
        else:
            locals()['ax'+str(i+1)].set_ylim(bottom=min(max(np.nanmin(obs_flux[M_obs]),vmin),np.nanmin(np.array(FL_conv[j])[M_FW]),locals()['ax'+str(i+1)].set_ylim()[0]),\
                                          top=max(min(np.nanmax(obs_flux[M_obs]),vmax,),np.nanmax(np.array(FL_conv[j])[M_FW]),locals()['ax'+str(i+1)].set_ylim()[1],1+0.2*ad))
        locals()['ax'+str(i+1)].set_xlim(xlim[i][0],xlim[i][1])
    for i in range(24):
        ad = np.diff(locals()['ax'+str(i+1)].set_ylim())
        l1,l2 = 1+0.1*ad,1+0.15*ad
        for t in range(len(transitions)):
            for cw in transitions[t]:
                if cw>=xlim[i][0] and cw<=xlim[i][1]:
                    locals()['ax'+str(i+1)].vlines(cw,l1,l2,colors=color[t])
                    
    custom_lines = [Line2D([0], [0], color=color[0], lw=1),
                    Line2D([0], [0], color=color[1], lw=1),
                    Line2D([0], [0], color=color[2], lw=1),
                    Line2D([0], [0], color=color[3], lw=1),
                    Line2D([0], [0], color=color[4], lw=1)]
    
    ax4.legend(custom_lines,['H','He','C','N','O'],bbox_to_anchor=(1.01, 1.05))  
    if models:
        prop_cycle = plt.rcParams['axes.prop_cycle']
        colors = prop_cycle.by_key()['color']
        MDL = []
        for l in range(len(models)):
            MDL.append(Line2D([0], [0], color=colors[l+1], lw=1))
        #ax24.legend(MDL,models,bbox_to_anchor=(1.01, 1.05))
        plt.figlegend(MDL,models,loc = (0.45,0.9))
    fig.supxlabel(r'$\lambda [\AA]$'),fig.suptitle(title)
    plt.show()
    return None

def fil_6(XY,title,color,transitions,vmax=1.1,vmin=0.6):
    XY_obs, XY_FW = XY
    wv,FL_conv,ot = XY_FW
    obs_wave,obs_flux,ot_1 = XY_obs
    
    fig = plt.figure(figsize=(12,18))
    ax1,ax5,ax10,ax14,ax18,ax22 = plt.subplot(6,4,1),plt.subplot(6,5,6), plt.subplot(6,4,9),plt.subplot(6,4,13),plt.subplot(6,4,17),plt.subplot(6,4,21)
    ax2,ax3,ax4 = plt.subplot(6,4,2,sharey =ax1),plt.subplot(6,4,3,sharey =ax1),plt.subplot(6,4,4,sharey =ax1)
    ax6,ax7,ax8,ax9=  plt.subplot(6,5,7,sharey =ax5),plt.subplot(6,5,8,sharey =ax5),plt.subplot(6,5,9,sharey =ax5),plt.subplot(6,5,10,sharey =ax5)
    ax11,ax12,ax13 = plt.subplot(6,4,10,sharey =ax10),plt.subplot(6,4,11,sharey =ax10),plt.subplot(6,4,12,sharey =ax10)
    ax15,ax16,ax17 = plt.subplot(6,4,14,sharey =ax14),plt.subplot(6,4,15,sharey =ax14),plt.subplot(6,4,16,sharey =ax14)
    ax19,ax20,ax21 = plt.subplot(6,4,18,sharey =ax18),plt.subplot(6,4,19,sharey =ax18),plt.subplot(6,4,20,sharey =ax18)
    ax23,ax24,ax25 = plt.subplot(6,4,22,sharey =ax22),plt.subplot(6,4,23,sharey =ax22),plt.subplot(6,4,24,sharey =ax22)
    xlim = [4079,4121],[4319,4361],[4839,4881],[6539,6581],\
        [4375,4405],[4910,4940],[4010,4040],[4455,4485],[4700,4730],\
        [4185,4215],[4525,4555],[4670,4700],[6670,6700],\
        [4622,4652],[4498,4528],[4043,4073],[4596,4626],\
        [4252,4282],[4632,4662],[5680,5710],[5786,5816],\
        [4060,4090],[4400,4430],[3946,3976],[5577,5607]
    lab = r'$H\delta+NIII4097$',r'$H\gamma$',r'$H\beta$',r'$H\alpha$',\
        r'$HeI4387$',r'$HeI4922$',r'$HeI+II4026$',r'$HeI4471$',r'$HeI4713$',\
        r'$HeI4200$',r'$HeII4541$',r'$HeII4686$',r'$HeI6678+HeII6683$',\
        r'$NII4630+NIII4634/40$',r'$NIII4511/15$',r'$NIV4058$',r'$NV4603/19$',\
        r'$CII4267$',r'$CIII4647$',r'$CIII5695$',r'$CIV5801$',\
        r'$OII4075$',r'$OII4414/16$',r'$OIII3961$',r'$OIII5592$'
    for i in range(25):
        M_FW = (wv>=xlim[i][0]) & (wv<=xlim[i][1])
        M_obs = (obs_wave>=xlim[i][0]) & (obs_wave<=xlim[i][1])
        if i==0 or i==4 or i==9 or i==13 or i==17 or i==21:
            locals()['ax'+str(i+1)].tick_params(which='both', direction='in',bottom=True,top=True,left=True,right=True)
        else:
            locals()['ax'+str(i+1)].tick_params(which='both', direction='in',bottom=True,top=True,left=True,labelleft=False,right=True)
        locals()['ax'+str(i+1)].xaxis.set_minor_locator(AutoMinorLocator(5))
        locals()['ax'+str(i+1)].yaxis.set_minor_locator(AutoMinorLocator(4))
        locals()['ax'+str(i+1)].plot(obs_wave[M_obs],obs_flux[M_obs])
        locals()['ax'+str(i+1)].plot(wv[M_FW],FL_conv[M_FW])
        locals()['ax'+str(i+1)].set_title(lab[i])
        ad = np.diff(locals()['ax'+str(i+1)].set_ylim())
        locals()['ax'+str(i+1)].set_ylim(bottom=min(max(np.nanmin(obs_flux[M_obs]),vmin),np.nanmin(FL_conv[M_FW]),locals()['ax'+str(i+1)].set_ylim()[0]),\
                                          top=max(1+0.25*ad,locals()['ax'+str(i+1)].set_ylim()[1],min(np.nanmax(FL_conv[M_FW]),vmax)))
        if i==0 or i==4 or i==9 or i==13 or i==16 or i==19 or i==22:
            locals()['ax'+str(i+1)].set_ylim(bottom=min(max(np.nanmin(obs_flux[M_obs]),vmin),np.nanmin(FL_conv[M_FW])),\
                                          top=max(min(np.nanmax(obs_flux[M_obs]),vmax,),np.nanmax(FL_conv[M_FW]),1+0.2*ad))            
        else:
            locals()['ax'+str(i+1)].set_ylim(bottom=min(max(np.nanmin(obs_flux[M_obs]),vmin),np.nanmin(FL_conv[M_FW]),locals()['ax'+str(i+1)].set_ylim()[0]),\
                                          top=max(min(np.nanmax(obs_flux[M_obs]),vmax,),np.nanmax(FL_conv[M_FW]),locals()['ax'+str(i+1)].set_ylim()[1],1+0.2*ad))
        locals()['ax'+str(i+1)].set_xlim(xlim[i][0],xlim[i][1])
    for i in range(25):
        ad = np.diff(locals()['ax'+str(i+1)].set_ylim())
        l1,l2 = 1+0.1*ad,1+0.15*ad
        for t in range(len(transitions)):
            for cw in transitions[t]:
                if cw>=xlim[i][0] and cw<=xlim[i][1]:
                    locals()['ax'+str(i+1)].vlines(cw,l1,l2,colors=color[t])
    
    custom_lines = [Line2D([0], [0], color=color[0], lw=1),
                    Line2D([0], [0], color=color[1], lw=1),
                    Line2D([0], [0], color=color[2], lw=1),
                    Line2D([0], [0], color=color[3], lw=1),
                    Line2D([0], [0], color=color[4], lw=1)]
    
    ax4.legend(custom_lines,['H','He','C','N','O'],bbox_to_anchor=(1.01, 1.05))  
    fig.supxlabel(r'$\lambda [\AA]$'),fig.suptitle(title)
    plt.show()
    return None


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

def spectrum_generation(model,elems,elems_abundance,VT,wv1=3800,wv2=7000,stp=0.01,\
                        folder_original_pkl=dir().modelgrid+dir().modelgrid.split('/')[-2]+'_pkl/original/'):
    '''
    Function to create the spectrum from a given element with the considered lines. 
    
    Parameters
    ----------
    model : str,
        Model name of the form 'T360g390He10Q130b10CNO' (INCLUDING the CNO specification, which is the only computed at this time)
    elems : str,
        Name of the elements to be added to the model (without CNO). At the moment: C, N, O, H_He
    elems_abundance : str
        Element abundance. It has to be in the pkl.
    VT : str,
        Microturbulence to be considered. It has to be in the pkl, and usually has the form 'VT007'.
    wv1 : float, optional
        First wavelength to be considered. The default is 3800.
    wv2 : float, optional
        Last wavelength to be considered. The default is 7000.
    stp : float, optional
        step size for the wavelength base. The default is 0.01.
    folder_original_pkl : str, optional
        Folder to find the pkl files from the different elements and models. The default is dir().modelgrid+dir().modelgrid.split('/')[-2]+'_pkl/original/'.

    Returns
    -------
    wave_range : 
        Wavelength array 
    flux_range
        Flux array (normalized to 1)

    '''
    wave_range = np.around(np.arange(wv1,wv2,stp),2)
    flux_range = np.zeros(len(wave_range))
    for elem in elems:
        with open(folder_original_pkl + model[:-3] + elem + '_lines.pkl', 'rb') as inp:
            profiles = pickle.load(inp, encoding='latin1')
            inp.close()
        LINES = pd.read_csv('../grid_components/lines_'+elem+'.dat', header=None,sep='\s+')
        metallicity = elems_abundance[elem]
        for line in LINES[0]:
            if line in considered_lines:
                line_wavelength = profiles.get_profile(line,str(metallicity),VT).wvl 
                if any(line_wavelength>=wv1) and any(line_wavelength<=wv2):      
                    line_flux = profiles.get_profile(line,str(metallicity),VT).flx 
                    print(line,max(line_flux),line_wavelength[np.argmax(line_flux)])
                    f1 = interp1d(line_wavelength,line_flux-1, kind='linear')  
                    mn, Mx = np.nanmin(line_wavelength),np.nanmax(line_wavelength)
                    mask = (wave_range>=mn) & (wave_range<=Mx)
                    masked_wave = wave_range[mask]
                    flx = f1(masked_wave)
                    flux_range[mask] += flx
    return wave_range,flux_range+1

 
def synthetic_model(star,model = None,beta = 1.5):
    '''
    Function to generate the complete synthetic spectra (broadened, etc)

    Parameters
    ----------
    star : str
        Name of the star 
    model : str, optional
        Name of the model to consider. The default is None.
        AT THE MOMENT, if model is None, it uses the one selected according to Carneiro, and will be implemented the one closest to Holgado
    beta : float, optional
        Darkening parameter in the rotational broadening. It has to be between 0 and 1.5. The default is 1.5.

    Returns
    -------
    TYPE
        DESCRIPTION.
    TYPE
        DESCRIPTION.

    '''
    E_C, E_N, E_O, vmic = Model[star][4:]
    E_H_He = E_N
    VT = find_nearest(FW_mic,vmic)
    VT = 'VT000'[:-len(str(int(np.rint(VT))))]+str(int(np.rint(VT)))
    elems = 'C','N','O','H_He'
    abundances = {}
    for elem in elems:
        metallicities = pd.read_csv('../grid_components/abundances_'+elem+'.dat', header=None,sep='\s+')[1]
        abundances[elem] = find_nearest(metallicities, locals()['E_'+elem]-12)
    
    # model_10, model_15 = best_models[star][:-4]+'0CNO',best_models[star][:-4]+'5CNO'
    # Now we have to compute the broadening
    if not model:
        model =  best_models[star][:-4]+'0CNO'
    vsini, vmac = broadening[star]
    #for model in [model_10, model_15]: At the moment we will fix beta = 1.0
    synthetic = spectrum_generation(model,elems,abundances,VT)
    FL_conv = np.zeros(len(synthetic[1]))
    broadening_sections = pd.read_csv('line_broadening_divisions.txt',skiprows=1, header=None,sep='\s+')
    
    for i in broadening_sections.index:
        m = broadening_sections[0][i],broadening_sections[1][i]
        l0 = np.mean(m)
        M = (synthetic[0]>m[0]) & (synthetic[0]<m[1])
        WAVE, FLUX = synthetic[0][M],synthetic[1][M]
        WV, FL = lb.final_resolution(WAVE, FLUX-1, l0, vsini, vmac, R, beta=beta)
        f1 = interp1d(WV, FL)
        X = np.arange(start=WAVE[0],stop=WAVE[-1]+np.diff(WAVE)[0]-1e-3,step=np.diff(WAVE)[0])
        Y = f1(X) 
        FL_conv[M] += Y
    FL_conv = FL_conv
    return synthetic[0],FL_conv+1
    