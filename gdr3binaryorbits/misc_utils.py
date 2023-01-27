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