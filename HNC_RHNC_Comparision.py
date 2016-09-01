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

T = 0.90
rho = 0.70 
Ur_ref = Potentials.WCAPotential(r,sig=1.0, eps=1.0, m=12, n=6)
Ur = Potentials.LJPotential(r,sig=1.0, eps=1.0, m=12, n=6)

Fig = plt.figure()
Fig1 = Fig.add_subplot(121)
Fig2 = Fig.add_subplot(122)
#HNC solver:
Flag, hr, cr, er, hk, Sk = HNC.OZSolver_HNC_Iterative(r, k, Ur, rho, kT=T, maxiter=10000, w_old_start=0.50,w_old_max=0.99,tol=1e-8, cr_guess=cr_guess)    
Fig1.plot(r,hr, 'b-',linewidth=3)
Fig2.plot(k,Sk, 'b-',linewidth=3)

#RHNC Solver
Flag, hr, cr, er, hk, Sk = RHNC.OZSolver_RHNC_Iterative(r, k, Ur, Ur_ref, rho, kT=T, maxiter=10000, w_old_start=0.50,w_old_max=0.99,tol=1e-8, cr_guess=cr_guess)    
Fig1.plot(r,hr, 'r-', linewidth=2)
Fig2.plot(k,Sk, 'r-', linewidth=2)

Fig1.set_xlim([0,3*sig])
#Fig2.set_xlim([0,30])
Fig2.set_xscale('log')

plt.show()