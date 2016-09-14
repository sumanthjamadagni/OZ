# -*- coding: utf-8 -*-

import numpy as np
from scipy.integrate import simps 
import sys 

def LJPotential_NonDim(r, T_LJ=1.0):
    '''
    Sigma is 1 - i.e., r is in reduced units of r/sigma
    '''
    U_LJ = 4 * ((1./r)**12 - (1./r)**6) / T_LJ
    return U_LJ



def SoftSpherePotential(r, sig=1.0, eps=1.0, m=12):
    U_SS = 4 * eps * ((sig/r)**m) 
    return U_SS



def LJPotential(r, sig=1.0, eps=1.0, m=12, n=6):
    if m < n:
        print "Wrong LJ type potential. m < n!!"
        print "Exiting"
        sys.exit()
    U_LJ = 4 * eps * ((sig/r)**m - (sig/r)**n) 
    return U_LJ


def WCAPotential(r, sig=1.0, eps=1.0, m=12, n=6):
    p = 1.0/(m-n)
    rcut = sig * (float(m)/float(n))**p
    U_temp = LJPotential(r, sig=sig, eps=eps) + eps
    U_WCA = np.where(r < rcut, U_temp, 0)
    return U_WCA


def SALRPotential(r, sig=1.0, eps=1.0, d=1.0, A=2.0, m=12, n=6):
    U_SALR = LJPotential(r,sig=sig,eps=eps,  m=m, n=n) + A*np.exp(-(r-sig)/d)/r
    return U_SALR

def SALRPotential_v2(r, sig=1.0, eps_rep=1.0, eps_att=1.0, d=1.0, A=2.0, m=12, n=6):
    U_SALR = 4*(eps_rep*(sig/r)**m - eps_att*(sig/r)**n)  + A*np.exp(-(r-sig)/d)/r
    return U_SALR

def SRLRPotential(r, sig=1.0, eps=1.0, d=1.0, A=2.0, m=12, n=6):       
    U_SRLR = WCAPotential(r,sig=sig,eps=eps, m=m, n=n) + A*np.exp(-(r-sig)/d)/r
    return U_SRLR


def DPDPotential(r, A=15.0):
    U_DPD = np.where(r < 1.0, 0.5*A*(1-r)**2, 0.0)
    return U_DPD


#Second virial coefficient
def CalcB2(r,Ur, kT=1.0):
    I = (np.exp(-Ur/kT)-1.0) * r**2
    B2 = -2*np.pi * simps(I,r)

    return B2

