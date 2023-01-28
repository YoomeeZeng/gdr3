# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import corner
from matplotlib import gridspec
from .misc_utils import *

def set_rcparams():
    plt.rcParams.update({'font.size': 10})
    plt.rcParams.update({'xtick.major.width': 1})
    plt.rcParams.update({'ytick.major.width': 1})
    plt.rcParams.update({'xtick.minor.width': 1})
    plt.rcParams.update({'ytick.minor.width': 1})     
    plt.rcParams.update({'xtick.major.size': 8})
    plt.rcParams.update({'ytick.major.size': 6})
    plt.rcParams.update({'xtick.minor.size': 3})
    plt.rcParams.update({'ytick.minor.size': 3})                 
    plt.rcParams.update({'font.family': 'STIXGeneral'})
    plt.rcParams.update({'mathtext.fontset': 'cm'})    

set_rcparams()

def plot_sb1_rv(rv_df,params_dict,rv_samples=pd.DataFrame()):
    
    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(20,10))
    
    for axis in ['top','bottom','left','right']:
                    ax.spines[axis].set_linewidth(1.5)     
                    
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    ax.set_xlabel(r'$\rm Orbital~Phase$',fontsize=54)  
    ax.set_ylabel(r'$\rm RV~ [km/s]$',fontsize=54)     
    ax.tick_params(direction='in',axis='both',which='minor',length=3,width=2,labelsize=36)
    ax.tick_params(direction='in',axis='both',which='major',length=6,width=2,labelsize=36)
    ax.minorticks_on()   
    

    if len(rv_samples)>0:
        
        for i in range(len(rv_samples)):
            rv, phase=get_sb1_orbit_samples(rv_samples.period[i],rv_samples.ecc[i],rv_samples.K1[i],rv_samples.arg_per[i],rv_samples.gamma[i],rv_samples.t_peri[i])
            
            sampled=pd.DataFrame({'phase':phase,'RV':rv})
            sampled.sort_values(by='phase',inplace=True)
            ax.plot(sampled.phase,sampled.RV,c='k',lw=0.5,alpha=0.1)
            ax.plot(sampled.phase+1,sampled.RV,c='k',lw=0.5,alpha=0.1)
        

    ax.plot(rv_df.phase,rv_df.RV,c='r',lw=3,alpha=1,label=r'$\rm Gaia~ DR3~ SB1~ Orbit$')
    ax.plot(rv_df.phase+1,rv_df.RV,c='r',lw=3,alpha=1)   
    
    gsrc=params_dict['gdr3_source']
    plt.title(f'Gaia DR3 {gsrc}',fontsize=48)
    
    ax.legend(fontsize=26,loc=1)
    ax.set_xlim(-0.05,2.05)
    gsrc=params_dict['gdr3_source']
    plt.title(f'Gaia DR3 {gsrc}',fontsize=48)
    plt.show()
    
def plot_sb1_orbit_data_comparison(rv_df,data,params_dict,rv_samples=pd.DataFrame()):
    
    fig, axes = plt.subplots(nrows=2, ncols=1,figsize=(25,25))
    gs = gridspec.GridSpec(2, 1, height_ratios=[2,1]) 
    
    ax1=plt.subplot(gs[1])
    ax=plt.subplot(gs[0], sharex=ax1)
    
    for axis in ['top','bottom','left','right']:
                    ax.spines[axis].set_linewidth(1.5) 
                    
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    ax.set_ylabel(r'$\rm RV~ [km/s]$',fontsize=54)     
    ax.tick_params(direction='in',axis='both',which='minor',length=3,width=2,labelsize=36)
    ax.tick_params(direction='in',axis='both',which='major',length=6,width=2,labelsize=36)
    ax.minorticks_on()   
    
    markers=['.','s','d']
    colors=['dodgerblue','C2','C3','C4']
    marker_sizes=[40,20,25]
    
    i=0
    
    #Plot the data
    for n,gr in data.groupby('obs_src'):
        
        ax.errorbar(gr.phase,gr.RV,yerr=gr.RV_err,markersize=marker_sizes[i],linestyle='none',label=n,marker=markers[i],alpha=0.75,c=colors[i])  
        ax.errorbar(gr.phase+1,gr.RV,yerr=gr.RV_err,markersize=marker_sizes[i],linestyle='none',marker=markers[i],alpha=0.75,c=colors[i]) 
            
        i+=1
    
    #Plot the orbit samples
    if len(rv_samples)>0:
        
        for i in range(len(rv_samples)):
            rv, phase=get_sb1_orbit_samples(rv_samples.period[i],rv_samples.ecc[i],rv_samples.K1[i],rv_samples.arg_per[i],rv_samples.gamma[i],rv_samples.t_peri[i])
            
            sampled=pd.DataFrame({'phase':phase,'RV':rv})
            sampled.sort_values(by='phase',inplace=True)
            ax.plot(sampled.phase,sampled.RV,c='k',lw=0.5,alpha=0.1)
            ax.plot(sampled.phase+1,sampled.RV,c='k',lw=0.5,alpha=0.1)
        
    #Plot the published model
    ax.plot(rv_df.phase,rv_df.RV,c='r',lw=3,alpha=1,label=r'$\rm Gaia~ DR3~ SB1~ Orbit$')
    ax.plot(rv_df.phase+1,rv_df.RV,c='r',lw=3,alpha=1)  
    ax.legend(fontsize=36)

    #Now plot residuals with the Gaia model    
    ax1.yaxis.set_ticks_position('both')
    ax1.xaxis.set_ticks_position('both')  
    ax1.tick_params(direction='in',axis='both',which='minor',length=3,width=2,labelsize=36)
    ax1.tick_params(direction='in',axis='both',which='major',length=6,width=2,labelsize=36)
    ax1.minorticks_on()  
    ax1.set_xlabel(r'$\rm Orbital~Phase$',fontsize=54)  
    
    for axis in ['top','bottom','left','right']:
        ax1.spines[axis].set_linewidth(1.5)  
    
    
    #Plot the residuals
    i=0
    plt.gca().set_prop_cycle(None)
    for n,gr in data.groupby('obs_src'):
        
        ax1.errorbar(gr.phase,gr.RV_resid,yerr=gr.RV_err,markersize=marker_sizes[i],linestyle='none',label=n,marker=markers[i],alpha=0.75,c=colors[i])  
        ax1.errorbar(gr.phase+1,gr.RV_resid,yerr=gr.RV_err,markersize=marker_sizes[i],linestyle='none',label=n,marker=markers[i],alpha=0.75,c=colors[i]) 
            
        i+=1
            
    ax1.set_ylabel(r'$\rm Residuals\,[km/s]$',fontsize=54)
    
    #ax1.set_ylim(-1.5,1.5)
    ax.set_xlim(-0.05,2.05)
    ax1.set_xlim(-0.05,2.05)
    ax1.axhline(0,c='k',alpha=0.5,ls='--',lw=3)
    
    #ax.text(0.9,-15,r'$P_{\rm orb}=81.16839^{+0.00044}_{-0.00045}~\rm d$',fontsize=42)
    gsrc=params_dict['gdr3_source']
    plt.title(f'Gaia DR3 {gsrc}',fontsize=48)   
    
    plt.show()
    