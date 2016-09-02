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

T = 0.80
rho = 0.75 
Ur_ref = Potentials.WCAPotential(r,sig=1.0, eps=1.0, m=12, n=6) #,A=1.0,d=2.0)
Ur = Potentials.SALRPotential(r,sig=1.0, eps=1.0, m=12, n=6, A=0.0,d=2.0)

Uk = OZF.FourierBesselTransform(Ur,r,k)
plt.plot(k,Uk,'r-d')

Fig = plt.figure(figsize=[15,12])
Fig1 = Fig.add_subplot(221)
Fig2 = Fig.add_subplot(222)
Fig3 = Fig.add_subplot(223)
Fig4 = Fig.add_subplot(224)

Fig1.set_xlabel('r', fontsize=FS)
Fig2.set_xlabel('k', fontsize=FS)
Fig3.set_xlabel('r', fontsize=FS)
Fig4.set_xlabel('r', fontsize=FS)

Fig1.set_ylabel('g(r)', fontsize=FS)
Fig2.set_ylabel('S(k)', fontsize=FS)
Fig3.set_ylabel('c(r)', fontsize=FS)
Fig4.set_ylabel('|rh(r)|', fontsize=FS)


cr_guess = None


#HNC solver:
Flag, hr_ref, cr_ref, er_ref, hk_ref, Sk_ref = HNC.OZSolver_HNC_Iterative(r, k, Ur_ref, rho, kT=T, maxiter=10000, w_old_start=0.50,w_old_max=0.99,tol=1e-8, cr_guess=cr_guess)    
Fig1.plot(r,hr_ref, 'g-',linewidth=4, label='$U_{rep}$, HNC')
Fig2.plot(k,Sk_ref, 'g-',linewidth=4)
Fig3.plot(r, cr_ref, 'g-', linewidth=4)
Fig4.plot(r, np.abs(hr_ref*r), 'g-', linewidth=4)
cr_guess = cr_ref

#HNC solver:
Flag, hr, cr, er, hk, Sk = HNC.OZSolver_HNC_Iterative(r, k, Ur, rho, kT=T, maxiter=10000, w_old_start=0.50,w_old_max=0.99,tol=1e-8, cr_guess=cr_guess)    
Fig1.plot(r,hr, 'b-',linewidth=3,label='$U_{att}$, HNC')
Fig2.plot(k,Sk, 'b-',linewidth=3)
Fig3.plot(r, cr, 'b-', linewidth=3)
Fig4.plot(r, np.abs(hr*r), 'b-', linewidth=3)

#RHNC Solver

Flag, hr, cr, er, hk, Sk = RHNC.OZSolver_RHNC_Iterative(r, k, Ur, Ur_ref, rho, kT=T, maxiter=10000, w_old_start=0.50,w_old_max=0.99,tol=1e-8, cr_guess=cr_guess)    
Fig1.plot(r,hr, 'r-', linewidth=2, label='$U_{att}$, RHNC')
Fig2.plot(k,Sk, 'r-', linewidth=2)
Fig3.plot(r, cr, 'r-', linewidth=2)
Fig4.plot(r, np.abs(hr*r), 'r-', linewidth=2)

Fig1.set_xlim([0,3*sig])
Fig2.set_xscale('log')
Fig3.set_xlim([0,2*sig])

Fig4.set_xscale('log')
Fig4.set_yscale('log')
Fig4.set_xlim([1.0, r[-1]])

Fig1.legend(loc='lower right')

plt.show()