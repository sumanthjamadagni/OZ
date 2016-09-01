# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.fftpack import dst, idst #discrete sine and inverse sine transforms
from collections import OrderedDict
from scipy.integrate import simps 

#------------------------
Na =  6.023e23 #/mol
R = 8.314e-3 #kJ/mol/K

pi = np.pi 
DIM=3 
#nr = 2048
#dr = 0.02
#dk = np.pi/nr/dr  #As Patrick Warren describes. 

#Normalization factor for the fourier transforms 
FAC_FT = 1.0 
FAC_INV_FT = 1.0/(2.0 * pi)**DIM 

#Normalization factor for fourier sine transfrom in 3D
FAC_FST = 4*pi
FAC_INV_FST = 4*pi

FAC_DFT = 1./2.#DON"T UNDERSTAND ORIGIN OF THIS EXACTLY. IT'S BECAUSE OF SINE TRANSFORM WHICH IS ONLY FOR ODD FUNCTIONS I THINK
FAC_INV_DFT = 1./2.#DON"T UNDERSTAND ORIGIN OF THIS EXACTLY. 

FAC_FBT = FAC_FT * FAC_FST * FAC_DFT
FAC_INV_FBT = FAC_INV_FT * FAC_INV_FST *FAC_INV_DFT

def Create_r_k(dr=0.02, nr=1024):
    r = np.arange(nr)*dr + dr #Can't have 0 in thr r array 
    dk = np.pi/nr/dr  #As Patrick Warren describes. 
    k = np.arange(nr)*dk + dk #Can't have 0 in the k array
    
    return r, k

def cr_initguess(Ur, kT=1.0, type='zeros'):
    if type == 'zeros':
        cr = np.zeros(len(Ur))
    elif type == 'MSA':
        cr = np.exp(-Ur/kT) - 1.0
    else:
        print "Unknown type for cr initial guess"
        print "options: zeros or MSA"
    return cr

def FourierBesselTransform(f, r, k):
    #Eqn. 2 of Patrick Warren's documentation
    dr = r[2]-r[1]
    fk = FAC_FBT * dst(f * r)/k * dr 
    return fk


def InvFourierBesselTransform(fk, r, k):
    #Eqn. 4 of Patrick Warren's documentation, using the inverse
    #discrete sine transform for the inv fourier bessel transform    
    dk = k[2]-k[1]    
    f =  FAC_INV_FBT * idst(fk * k)/r * dk         
    return f

def Calculate_ek(ck, rho):
    #Eqn 8 in Patrick Warren's documentation
    ek = ck/(1-rho*ck) - ck
    return ek


def GPY_Closure(gr, Ur, kT=1.0):
    #Generalized Percus Yevick closure for continuous potentials. 
    cr = gr * (1.0  - np.exp(Ur/kT))       
    return cr

def RPA_Closure(Ur, kT=1.0):
    #RPA closure for soft potentials. 
    #Conditions for mathematical convergence (for ck to exist): 
    #  1.  Must decay faster than 1/r^2 for large r 
    #  2. diverge slower than 1/r^2 for small r.  
    cr = -Ur/kT      
    return cr

def Calculate_hr(er, cr):
    hr = er + cr
    return hr

def Calculate_gr(hr):
    gr = hr + 1
    return gr 
    
def ListToTabbedStr(List):
    tab = "\t"
    N = len(List)
    strval = str(List[0]) + tab
    for i in range(1,N-1):
        strval = strval + str(List[i]) + tab 
    
    strval = strval + str(List[N-1]) + "\n" 
    return strval

def PY_HS_Analytical(r, rho, sig=1.0):
    '''
    From: http://www.sklogwiki.org/SklogWiki/index.php/Exact_solution_of_the_Percus_Yevick_integral_equation_for_hard_spheres
    '''
    dr = r[2] - r[1]
    nr = len(r)
    k = np.pi*np.arange(1,nr+1)/(dr*nr)
    
    eta = np.pi/6.0 * sig**3 * rho 
    
    denom = (1.0 - eta)**4
    cr = (1. + 2*eta)**2 - 6*eta*(1.0 + 0.5*eta)**2*(r/sig) + eta*(1.+2.*eta)**2 * (r/sig)**3/2.0
    cr = np.where(r < sig, -cr/denom, 0)
    
    ck = FourierBesselTransform(cr,r,k)
    hk = ck / (1.0 - rho * ck)
    hr = InvFourierBesselTransform(hk,r,k)
    hr = np.where(r < sig, -1, hr)
    
    return cr, hr

