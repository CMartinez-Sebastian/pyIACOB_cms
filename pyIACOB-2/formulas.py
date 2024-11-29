#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 27 11:06:23 2021

@author: charlie
"""

import numpy as np 
import matplotlib.pyplot as plt
from astropy.io import fits
import math as mth
import os
from scipy.optimize import curve_fit
from scipy import integrate
plt.close('all')
plt.rcParams['figure.figsize'] = [8, 6]
c = 299792.458
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx],idx
def ajuste_mas(fit):
    hdul = fits.open(fit)
    #y = fits.getdata(fit)[0]
    head=hdul[0].header.copy()
    y=hdul[0].data.copy()
    first = head['CRVAL1']
    step = head['CDELT1']
    vel = head['I-VBAR']
    beta = vel/c# 1 LEER ESPECTRO, A CONTINUACIÓN CORRECCIÓN 
    x2= np.linspace(first,first+len(y[0])*step,len(y[0]))
    x2 = x2/(1+beta)
    hdul.close()
    return x2, y[0]
def ajuste_menos(fit):
    hdul = fits.open(fit)
    #y = fits.getdata(fit)[0]
    head=hdul[0].header.copy()
    y=hdul[0].data.copy()
    first = head['CRVAL1']
    step = head['CDELT1']
    vel = head['I-VBAR']
    beta = vel/c
    x2= np.linspace(first,first+len(y[0])*step,len(y[0]))
    x2 = x2/(1-beta)
    hdul.close()
    return x2, y[0]
def gaussian_flux(x, a, sigma, mu):
    # Function to create a flux as a gaussian. It it normalized at 1. 
    r = 1+a/(sigma*np.sqrt(2*np.pi))*np.exp(-(x-mu)**2/(2*sigma**2))
    return(r)  
def adjust(function, data, pars0=None):
    # Function to adjust 3 values (a, sigma, mu). To do a more clear code. 
    pars, cov = curve_fit(f=function, xdata=data[0], ydata=data[1], p0=pars0, bounds=(-np.inf, np.inf))
    return(pars, cov) 
def Pars0(dat):   
    # Let's calculate here a first approximation for the initial values to iniciate 
    # the fit. It could improve the computing time.
    Mu = sum(dat[0])/len(dat[0])
    sigma2 = sum(dat[0]**2) /len(dat[0]) - Mu**2
    Sig = sigma2**0.5   
    Max = max(dat[1]-1)*Sig*np.sqrt(2*np.pi)    
    pars0 = [Max, Sig, Mu]
    return(pars0)
def fit_gaussian(central_ojo,ristra):
    x,y=ristra[0],ristra[1]
    E=find_nearest(x,central_ojo)
    E_min=find_nearest(y,min((y[E[1]-100:E[1]+100])))
    E_max=find_nearest(y,max((y[E[1]-100:E[1]+100])))
    if (E_max[0]-1)>(1-E_min[0]):
        E_min=E_max
    mask=abs(y[E_min[1]-400:E_min[1]+400]-1)>=0.15*(abs(1-E_min[0]))
    y_E=y[E_min[1]-400:E_min[1]+400]*mask
    x_E=x[E_min[1]-400:E_min[1]+400]*mask
    
    # Ajuste de gaussiana al perfil 
    X_E=[x_E[x_E !=0],y_E[y_E !=0]]
    pairs, covar=adjust(gaussian_flux,X_E,Pars0(X_E))
    x_lin=np.linspace(x_E[x_E !=0][0],x_E[x_E !=0][-1],100)
    return [x_lin,gaussian_flux(x_lin,*pairs),pairs]

def elem(fich,table):
    fichero= open(fich,'r')
    for lin in fichero:
        l = lin.split()
        try:
            table.append(float(l[0]))
        except:
            np.nan
    fichero.close()
    return table
def get_info(fit,key):
    hdul = fits.open(fit)
    #y = fits.getdata(fit)[0]
    head=hdul[0].header.copy()
    info=head[key]
    del head
    hdul.close()
    return info
def minimo(central_ojo,ristra):
    x,y=ristra[0],ristra[1]
    E=find_nearest(x,central_ojo)
    A = max(abs(1-y[E[1]-100:E[1]+100]))
    E_min=find_nearest(y,min((y[E[1]-100:E[1]+100])))
    E_max=find_nearest(y,max((y[E[1]-100:E[1]+100])))
    if (E_max[0]-1)>(1-E_min[0]):
        E_min=E_max
    mask=abs(y[E_min[1]-400:E_min[1]+400]-1)>=0.15*(abs(1-E_min[0]))#y[E_min[1]-400:E_min[1]+400]<=0.93
    return E_min

def zero_moment(x_l,y_l,line):
    v_l = (x_l-line)*c/line
    F_l=1-y_l
    y_int0=integrate.cumtrapz(F_l, v_l, initial=0)    
    return y_int0[-1]    

def first_moment(x_l,y_l,line):
    #globals()[str(BJD)+'_BJDSiIII'],globals()[str(BJD)+'_stdSiIII'],globals()[str(BJD)+'_xSiIII'],globals()[str(BJD)+'_ySiIII']=BJD,std,x_l,y_l
    #BJD_list.append(BJD)
    v_l = (x_l-line)*c/line
    F_l=1-y_l
    y_int1 = integrate.cumtrapz(F_l*v_l, v_l, initial=0)
    return y_int1[-1]

def radial_velocity(x_l,y_l,line):
    v_l = (x_l-line)*c/line    
    return fit_gaussian(0,[v_l,y_l])[-1][-1]

def second_moment(x_l,y_l,line,centre=None):
    #globals()[str(BJD)+'_BJDSiIII'],globals()[str(BJD)+'_stdSiIII'],globals()[str(BJD)+'_xSiIII'],globals()[str(BJD)+'_ySiIII']=BJD,std,x_l,y_l
    #BJD_list.append(BJD)
    v_l = (x_l-line)*c/line-centre
    F_l=1-y_l
    y_int2 = integrate.cumtrapz(F_l*v_l**2, v_l, initial=0)
    return y_int2[-1]
def third_moment(x_l,y_l,line,centre=None):
    #globals()[str(BJD)+'_BJDSiIII'],globals()[str(BJD)+'_stdSiIII'],globals()[str(BJD)+'_xSiIII'],globals()[str(BJD)+'_ySiIII']=BJD,std,x_l,y_l
    #BJD_list.append(BJD)
    v_l = (x_l-line)*c/line-centre
    F_l=1-y_l
    y_int3 = integrate.cumtrapz(F_l*v_l**3, v_l, initial=0)
    return y_int3[-1]

def PCA(dat,t,Vmin,Vmax,lam,tit):
    # t== Último autovector que considero. 
    V=dat.copy()
    D=np.reshape(V,(len(V[:,0,0])*len(V[0,:,0]),96))
    # Saco la matriz de correlación
    C=np.matrix(D.transpose())*np.matrix(D)
    #Diagonalizamos y ordenamos 
    evals, evecs = np.linalg.eig(C)
    vector = []
    for i in range(0,len(evecs)):
        vector.append(evecs[:,i])
    Z=sorted(zip(evals, vector),reverse=True)
    evals,evecs=[Z[i][0] for i in range(0,len(evals))],[Z[i][1] for i in range(0,len(evals))]
    # Sacamos la matriz de coeficientes en la nueva base, diagonalizada 
    B=evecs.copy()
    # Sacamos la matriz de coeficientes en la nueva base, diagonalizada 
    A=np.dot(D,B)[:,:,0]
    D_PCA=np.dot(A[:,:t],np.squeeze(B)[:t])
    # Redimensionamos la matriz de datos, así como la de residuos. 
    D_PCA_res=D-D_PCA    
    res_V1=np.reshape(D_PCA,(len(V[:,0,0]),len(V[0,:,0]),96))
    res_V2=np.reshape(D_PCA_res,(len(V[:,0,0]),len(V[0,:,0]),96))
    # Ponemos el límite del colormap, de tal modo que TODAS las gráficas tengan la 
    # misma escala y se facilite el estudio de los residuos. 
    # Pintamos las 3 imágenes con la misma escala en una misma figura para comparar. 
    '''
    fig, axes = plt.subplots(1, 3)
    ax1, ax2,ax3=axes[0],axes[1],axes[2]
    fig.suptitle('data cleansing for '+str(tit)+r', $\lambda$='+str(lamb_s[25])[:8]+r'$\AA$',fontsize=16)
    im1 = ax1.imshow(dat[:,:,lam],cmap='binary_r',vmax=Vmax,vmin=Vmin,extent=[0,len(V[0,:,0])*0.080,0,len(V[:,0,0])*0.080])
    ax2.imshow(res_V1[:,:,lam],cmap='binary_r',vmax=Vmax,vmin=Vmin,extent=[0,len(V[0,:,0])*0.080,0,len(V[:,0,0])*0.080])
    ax3.imshow(res_V2[:,:,lam],cmap='binary_r',vmax=Vmax,vmin=Vmin,extent=[0,len(V[0,:,0])*0.080,0,len(V[:,0,0])*0.080])
    ax1.set_title('Original',fontsize=18)
    ax2.set_title('Clean',fontsize=18)
    ax3.set_title('Residues',fontsize=18)
    ax1.set_xlabel('Distance from disk center [arcsec]',fontsize=18)
    ax1.set_ylabel('Distance from disk center [arcsec]',fontsize=18)
    ax2.set_xlabel('Distance from disk center [arcsec]',fontsize=18)
    ax2.set_ylabel('Distance from disk center [arcsec]',fontsize=18)
    ax3.set_xlabel('Distance from disk center [arcsec]',fontsize=18)
    ax3.set_ylabel('Distance from disk center [arcsec]',fontsize=18)
    ax1.tick_params(axis='y', labelsize=18)
    ax2.tick_params(axis='y', labelsize=18)
    ax3.tick_params(axis='y', labelsize=18)
    ax1.tick_params(axis='x', labelsize=18)
    ax2.tick_params(axis='x', labelsize=18)
    ax3.tick_params(axis='x', labelsize=18)
    color=fig.colorbar(im1,ax=axes.ravel().tolist(), shrink=0.55)
    color.ax.tick_params(labelsize=18)
    color.set_label('Circular polarization [V/I]', rotation=90,fontsize=18)
    plt.close('all')
    '''
    return(res_V1)    
    
    
    