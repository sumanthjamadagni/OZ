# -*- coding: utf-8 -*-
# Functions for post-processing results after OZ equations have been solved

import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.fftpack import dst, idst #discrete sine and inverse sine transforms
from collections import OrderedDict
from scipy.integrate import simps 

pi = np.pi

# --------------------------------------------
           
def Calc_Pvirial(r,gr,Ur, rho, kT=1.0, rmin=0.10):
    from scipy.integrate import simps 

    dr = r[2]-r[1]
    nmin = int(rmin/dr)
    dUdr = np.diff(Ur)/dr #of length len(r)-1 
    
    Integrand = 4*pi*r[:-1]**3 * gr[:-1]* dUdr
    Integral = simps(Integrand[nmin:],r[nmin:-1])
    Pv = rho * kT - rho**2/6.0  * Integral
    return Pv
    
def Calc_Pcompressibility(rho, kappa, opt='cum', offset=None, kT=1.0, B2=0.0):
    '''
    Both rho and kappa are 1D arrays.
    opt = cum
        Returns pressure as a function of rho - array containing N-1 elements. 
    opt = final
        Returns a scalar. 
   '''
    from scipy.integrate import cumtrapz, simps 
    
    if offset == None:
        try:
            P0 = (rho[0] + B2*rho[0]**2)*kT
        except IndexError:
            P0 = 0.0
    
    if opt == 'cum':
        P_compressibility = cumtrapz(1.0/(rho * kappa), rho) + P0
    else:
        P_compressibility = simps(1.0/rho*kappa, rho)
   
    return P_compressibility 
    
def CalculateKappa(Sk, rho, kT=1.0):
    kappa = Sk[0]/ rho / kT
    return kappa 
#--------------------------------------------------------------

def s2(r,gr,rho):
    '''
    2-particle excess entropy : Eqn. 9 in A Baranyi and DJ Evans, Phys. Rev. A., 1989 
    '''

    Integrand = np.where(gr > 0, -0.5 * rho * (gr * np.log(gr) - gr + 1.0), -0.5 * rho)
    s2 = rho * simps(Integrand,r)
    return s2

# --------------------------------------------------------------    
def CreateOutput(r, hr, cr, er, k, hk, Sk, OutFile):
    '''
    For diagnostics on a particular OZ solve run - write out total, direct and indirect correlations. 
    '''
    nr = len(r)
    OutFileH = open(OutFile, 'w')
    OutFileH.write('#r  h(r)  c(r)  e(r) k hk(k) Sk(k) \n')
    for i in range(nr):
        List = [r[i], hr[i], cr[i], er[i], k[i], hk[i], Sk[i] ]
        StrVal = ListToTabbedStr(List)
        OutFileH.write(StrVal)

#-----------------------------------------------
def CoordinationNumber(r,gr,rho, rcut=None):
    from scipy.integrate import cumtrapz
    #1 for the central particle itself. 
    nr = 1 + 4 * np.pi * rho * cumtrapz(r**2 * gr, r)
    if rcut:
        dr=r[2]-r[1]
        n = int(rcut/dr) + 1
        ncut = 1 + 4 * np.pi * rho * simps(r[:n]**2 * gr[:n], r[:n])
    if rcut:
        return nr, ncut
    else:
        return nr
    