# -*- coding: utf-8 -*-

# -*- coding: utf-8 -*-

# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt 
import matplotlib
from itertools import cycle

import Potentials
import OZ_Functions as OZF
import PP_Functions 
import FigFuncs
import HNC
import RHNC

colors = ['red', 'blue', 'green', 'cyan', 'black', 'orange']
ColorCycler = cycle(colors)

FS = 18
matplotlib.rc('xtick', labelsize=FS)
matplotlib.rc('ytick', labelsize=FS)

dr = 0.02
nr = 2048 
r,k = OZF.Create_r_k(dr, nr)
sig=1.0

T = 0.75
Ur_ref = Potentials.WCAPotential(r,sig=1.0, eps=1.0, m=12, n=6)
Ur = Potentials.LJPotential(r,sig=1.0, eps=1.0, m=12, n=6)

Fig = plt.figure(figsize=[12,10])
Fig1 = Fig.add_subplot(221)
Fig2 = Fig.add_subplot(222)
Fig3 = Fig.add_subplot(223)
Fig4 = Fig.add_subplot(224)

Fig1.set_xlabel('r', fontsize=FS)
Fig2.set_xlabel('k$\sigma$/2$\pi$', fontsize=FS)

Fig1.set_ylabel('g(r)', fontsize=FS)
Fig2.set_ylabel('S(k)', fontsize=FS)
Fig1.set_xlim([0,4*sig])
#Fig2.set_xlim([0,5])
Fig2.set_xscale('log')

B2 = Potentials.CalcB2(r,Ur,kT=T)

#Plotting U(r)
titleStr = "B2 = " + str(np.round(B2,3))
FigUr = FigFuncs.CreateFig(r,Ur, title=titleStr,xlabel='r', ylabel='U(r)', xleft=0, xright=4*sig, ybottom=np.min(Ur)-0.25, ytop=5.0, lw=2, marker='None')
try:
    FigUr_ref = FigFuncs.CreateFig(r,Ur_ref, title=titleStr,xlabel='r', ylabel='U(r)', xleft=0, xright=4*sig, ybottom=np.min(Ur)-0.25, ytop=5.0, lw=2, marker='None')
except NameError:
    pass
       
FigNR = plt.figure()
FigNR1= FigNR.add_subplot(111)
FigNR1.set_xlabel('r', fontsize=FS)
FigNR1.set_xlim([0.5,4])
FigNR1.set_yscale('log')
FigNR1.set_ylabel('n(r)', fontsize=FS)

#Plotting 2-particle excess entropy calculated.
FigS2 = plt.figure()
Figs2= FigS2.add_subplot(111)
Figs2.set_xlabel('density, $\\rho\sigma^3$', fontsize=FS)
Figs2.set_ylabel('s$_2$', fontsize=FS)


print "B2 = ", B2 

drho_vap = 0.005
rho_vap = np.arange(0.005, 0.10, drho_vap)

drho_liq = 0.01
rho_liq = np.arange(0.60, 0.90, drho_liq)
PlotFreq = 1

cr_guess = None
nrho_vap = len(rho_vap)
rho_converged_Vap  = []
Pressures_Vap = []
Compressibility_Vap = []
ExcessEntropy_Vap = []

nrho_liq = len(rho_liq)
rho_converged_Liq  = []
Pressures_Liq = []
Compressibility_Liq = []
ExcessEntropy_Liq = []

#rho = np.append(rho_vap)

cr_guess = None
for irho in range(nrho_liq):
    rho = rho_liq[irho]    
    Flag, hr, cr, er, hk, Sk = RHNC.OZSolver_RHNC_Iterative(r, k, Ur, Ur_ref, rho, kT=T, maxiter=10000, w_old_start=0.50,w_old_max=0.99,tol=1e-8, cr_guess=cr_guess)    
    if Flag == 0:
        cr_guess = cr #for next density
        P = PP_Functions.Calc_Pvirial(r,hr+1,Ur,rho,kT=T, rmin=0.25)
        print " P = ", P 
        Pressures_Liq.append(P)
        rho_converged_Liq.append(rho)
        kappa = PP_Functions.CalculateKappa(1+rho*hk, rho, kT=T)
        Compressibility_Liq.append(kappa)
        
        nr = PP_Functions.CoordinationNumber(r,hr+1,rho)
        s2 = PP_Functions.s2(r,hr+1,rho)
        ExcessEntropy_Liq.append(s2)
                        
        if irho % PlotFreq == 0:
            Fig1.plot(r, hr+1, linewidth=2, label='$\\rho = $' + str(rho))
            Fig2.plot(k*sig/(2*np.pi), 1 + rho*hk, linewidth=2)
            FigNR1.plot(r[1:], nr)
    else:
        pass
        


rho_converged_Liq = np.array(rho_converged_Liq)
P_gas = rho_converged_Liq * T + rho_converged_Liq**2 * T * B2 

P_compressibility = PP_Functions.Calc_Pcompressibility(rho_converged_Liq, Compressibility_Liq, opt='cum',kT=T, B2=B2)

Fig1.legend(loc='upper right', ncol=2, framealpha=0.50)
Fig3.set_xlabel('Density, $\\rho\sigma^3$', fontsize=FS)
Fig4.set_xlabel('Density, $\\rho\sigma^3$', fontsize=FS)
Fig3.set_ylabel('Pressure', fontsize=FS)
Fig4.set_ylabel('Compressibility', fontsize=FS)
Fig4.set_yscale('log')

Fig3.plot(rho_converged_Liq, Pressures_Liq, 'b-s', label='Virial')
Fig3.text(0.10, 0.70, "B$_2$=" + str(np.round(B2,2)), transform=Fig3.transAxes, fontsize=FS)
Fig3.plot(rho_converged_Liq, P_gas, 'g-d', label='IG + B2')
Fig3.plot(rho_converged_Liq[1:], P_compressibility, 'r-d', label='Compressibility')
Fig3.legend(loc='upper left')

Fig4.plot(rho_converged_Liq, Compressibility_Liq, 'r-d')
Fig.tight_layout()

Figs2.plot(rho_converged_Liq,ExcessEntropy_Liq, 'r-d')
