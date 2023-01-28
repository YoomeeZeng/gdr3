# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

def nr(g, dg, E, M,e,eps=1e-5):
    for i in range(50):
        E = E - (g(E,e,M))/dg(E,e,M)
    return E

def g(E,e,M):
    y=E-e*np.sin(E)-M
    return y

def dg(E,e,M):
    y=1-e*np.cos(E)
    return y


def rv(x,K,e,w,gamma):
    M=2*np.pi*(x)
    w=np.pi*w/180
    
    #Mean anomaly
    E=nr(g, dg, M,M,e,eps=1e-5)
    xx=np.sqrt((1+e)/(1-e))*np.tan(E/2)
    
    #nu is true anomaly
    nu=2*np.arctan(xx)
    y=K*(np.cos(nu+w)+e*np.cos(w))+gamma

    return y

def get_phased_sb1_rvs(phases,K1,ecc,arg_per,gamma):
    
    rvs=rv(phases,K1,ecc,arg_per,gamma)
    
    return rvs

def get_phased_sb2_rvs(phases,K1,ecc,arg_per,gamma,q):
    
    rvs_1=rv(phases,K1,ecc,arg_per,gamma)
    rvs_2=-rvs_1*q
    
    return rvs_1, rvs_2

def get_sb1_orbit_samples(period,ecc,K1,arg_per,gamma,t_peri):
    
    epochs=np.linspace(2457389.0-1000,2457389.0+1000,10000)
    phases=((epochs-t_peri)/period)%1  
    
    rv_sampled=get_phased_sb1_rvs(phases,K1,ecc,arg_per,gamma)
    
    return rv_sampled, phases
    

def get_phase(time,period,epoch):
    
    return ((time-epoch)/period)%1     

def get_predicted_rv(time,rv_samples):
    
    phases=np.vectorize(get_phase)(time,rv_samples.period,rv_samples.t_peri)
    rvs=np.vectorize(rv)(phases,rv_samples.K1,rv_samples.ecc,rv_samples.arg_per,rv_samples.gamma)
    
    rvs_median=np.median(rvs)
    rvs_std=np.std(rvs)
    
    return rvs_median, rvs_std
    
    
    
    