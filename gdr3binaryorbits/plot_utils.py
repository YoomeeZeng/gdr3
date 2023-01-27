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

def plot_rvs(rv_df):
    
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
    
    ax.plot(rv_df.phase,rv_df.RV,c='k',lw=3,alpha=0.5,label=r'$\rm Model$')
    ax.plot(rv_df.phase+1,rv_df.RV,c='k',lw=3,alpha=0.5)    
    ax.legend(fontsize=20,loc=2)
    ax.set_xlim(-0.05,2.05)
    plt.show()