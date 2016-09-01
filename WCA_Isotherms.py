# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt 
import Potentials
import OZ_Functions
import matplotlib
from itertools import cycle

colors = ['red', 'blue', 'green', 'cyan', 'black', 'orange']
ColorCycler = cycle(colors)

FS = 18
matplotlib.rc('xtick', labelsize=FS)
matplotlib.rc('ytick', labelsize=FS)

dr = 0.02
nr = 2048
r,k = OZ_Functions.Create_r_k(dr, nr)
Ur = Potentials.WCAPotential(r,sig=1.0, eps=1.0)


Fig = plt.figure()
Fig1 = Fig.add_subplot(221)
Fig2 = Fig.add_subplot(222)
Fig3 = Fig.add_subplot(223)
Fig4 = Fig.add_subplot(224)

Fig1.set_xlabel('r', fontsize=FS)
Fig2.set_xlabel('k', fontsize=FS)

Fig1.set_ylabel('g(r)', fontsize=FS)
Fig2.set_ylabel('S(k)', fontsize=FS)
Fig1.set_xlim([0,4])
Fig2.set_xlim([0,30])
T = 1.0
B2 = Potentials.CalcB2(r,Ur,kT=T)

nrho = 12
rhoArray = np.linspace(0.01,1.01,nrho)
cr_guess = None
rho_converged  = []
Pressures = []
Compressibility = []
for irho in range(nrho):
    rho = rhoArray[irho]
    Flag, hr, cr, er, hk, Sk = OZ_Functions.OZSolver_HNC_Iterative(r, k, Ur, rho, kT=T, maxiter=10000, w_old_start=0.50,w_old_max=0.99,tol=1e-8, cr_guess=cr_guess)
    cr_guess = cr
    if Flag == 0:
        P = OZ_Functions.Calc_Pvirial(r,hr+1,Ur,rho,kT=T, rmin=0.15)
        print " P = ", P 
        Pressures.append(P)
        rho_converged.append(rho)
        kappa = OZ_Functions.CalculateKappa(1+rho*hk, rho, kT=T)
        Compressibility.append(kappa)
        
        Fig1.plot(r, hr + 1.0)
        Fig2.plot(k, 1 + rho*hk)

rho_converged = np.array(rho_converged)
P_gas = rho_converged * T + rho_converged**2 * T * B2 

Fig3.set_xlabel('Density', fontsize=FS)
Fig4.set_xlabel('Density', fontsize=FS)
Fig3.set_ylabel('Pressure', fontsize=FS)
Fig4.set_ylabel('Compressibility', fontsize=FS)
Fig4.set_yscale('log')

Fig3.plot(rho_converged, Pressures, 'b-s')
Fig3.plot(rho_converged, P_gas, 'g-d')

Fig4.plot(rho_converged, Compressibility, 'r-d')

Fig.suptitle("T = " + str(T), fontsize=FS + 2)
Fig.tight_layout()