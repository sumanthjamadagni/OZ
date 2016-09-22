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
    U_temp = LJPotential(r, sig=sig, eps=eps, m=m, n=n)
    U_temp = U_temp - np.min(U_temp)
    U_WCA = np.where(r < rcut, U_temp, 0)
    return U_WCA


def LJPotential_v2(r, sig=1.0, eps=1.0, m=12, n=6):
    if m < n:
        print "Wrong LJ type potential. m < n!!"
        print "Exiting"
        sys.exit()
    U_LJ = 4 * ((sig/r)**m - eps*(sig/r)**n) 
    return U_LJ


def WCAPotential_v2(r, sig=1.0, eps=1.0, m=12, n=6):
    p = 1.0/(m-n)
    rcut = sig * (float(m)/float(n)/eps)**p
    U_temp = LJPotential(r, sig=sig, eps=eps, m=m, n=n)
    U_temp = U_temp - np.min(U_temp)
    U_WCA = np.where(r < rcut, U_temp, 0)
    return U_WCA


def SALRPotential(r, sig=1.0, eps=1.0, d=1.0, A=2.0, m=12, n=6):
    U_SALR = LJPotential(r,sig=sig,eps=eps,  m=m, n=n) + A*np.exp(-(r-sig)/d)/r
    return U_SALR



def SRLRPotential(r, sig=1.0, eps=1.0, d=1.0, A=2.0, m=12, n=6):       
    U_SRLR = WCAPotential(r,sig=sig,eps=eps, m=m, n=n) + A*np.exp(-(r-sig)/d)/r
    return U_SRLR

#---
def SALRPotential_v2(r, sig=1.0, eps=1.0, d=1.0, A=2.0, m=12, n=6):
    U_SALR = 4 * (sig/r)**m - 4 * eps*(sig/r)**n  + A*np.exp(-(r-sig)/d)/r
    return U_SALR

def SRLRPotential_v2(r, sig=1.0, eps=1.0, d=1.0, A=2.0, m=12, n=6):       
    U_SRLR = WCAPotential_v2(r,sig=sig,eps=eps, m=m, n=n) + A*np.exp(-(r-sig)/d)/r
    return U_SRLR

#---
def DPDPotential(r, A=15.0):
    U_DPD = np.where(r < 1.0, 0.5*A*(1-r)**2, 0.0)
    return U_DPD


#Second virial coefficient
def CalcB2(r,Ur, kT=1.0):
    I = (np.exp(-Ur/kT)-1.0) * r**2
    B2 = -2*np.pi * simps(I,r)

    return B2

#From Noro and Frenkel, JCP 2000: Extended Corresponding-State
#Behavior For Particles With Variable Range Attractions.

def Sigma_Eff(r,Urep, kT=1.0):
    #Eqn 4: Noro and Frenkel, JCP 2000
    I = 1.0 - np.exp(-Urep/kT)
    sigma = simps(I,r)
    return sigma


def CalcB2_reduced(r, Ur, Urep, kT=1.0):
    #Eqnation 6, Noro and Frenkel, JCP 2000
    B2 = CalcB2(r,Ur,kT)
    sigma_eff = Sigma_Eff(r,Urep, kT)

    B2_star = B2/(2.0/3.0 *np.pi * sigma_eff**3.0)
    return B2_star

def StickinessParameter(B2_star):
    #Eqnation 8, Noro and Frenkel, JCP 2000
    tau = 0.25/(1.0 - B2_star)
    return tau


def SWMapping(r, Ur, Urep, kT=1.0):

    eps = np.min(Ur)
    Tstar = -kT/eps #reduced temperature
    if Tstar < 0.0:
        print "Tstar is negative. No minima in potential!!"
        
    B2_star = CalcB2_reduced(r,Ur,Urep,kT)

    #Eqnation 10, Noro and Frenkel, JCP 2000
    lb = ((B2_star - 1.0)/(1.0 - np.exp(1.0/Tstar)) + 1.0)**(1.0/3.0)

    v = (B2_star - 1.0)/(1.0 - np.exp(1.0/Tstar)) + 1.0

    tau = 1.0 / ( 4.0 * (lb**3 - 1) * (np.exp(1.0/Tstar)-1) )

    return lb, tau , Tstar
    
    

    

    
    
    
    





