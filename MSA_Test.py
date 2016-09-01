# -*- coding: utf-8 -*-


#Comparing g(r) and S(k) for hard spheres - 
#from numerically solving the MSA closure and 
#comparing with the analutical Percus-Yevick result. 

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

dr = 0.021
nr = 2048
r,k = OZ_Functions.Create_r_k(dr, nr)
Ur = np.zeros(nr) #Hard spheres. 
kT  = 1.0
sig = 1.0

rho_list = [0.55, 0.58, 0.62] #[0.10, 0.20, 0.30, 0.40, 0.45, 0.48, 0.50, 0.55, 0.58, 0.62] #already in reduced units, rho*sig**3. 
rho_list = np.array(rho_list) 
nrho = len(rho_list)
eta_list = np.pi/6.0 * rho_list #packing fraction

Fig = plt.figure(figsize=[14,5])
Fig1 = Fig.add_subplot(131)
Fig2 = Fig.add_subplot(132)
Fig3 = Fig.add_subplot(133)

Fig1.set_xlim([0.5, 3])
Fig2.set_xlim([0, 1.5])
Fig3.set_xlim([0,20])

Fig1.set_title('gr')
Fig2.set_title('cr')
Fig3.set_title('Sk')

for irho in range(nrho):
    lc = next(ColorCycler)
    rho = rho_list[irho]; eta = eta_list[irho]
    
    cr, hr = OZ_Functions.PY_HS_Analytical(r, rho)
    hk = OZ_Functions.FourierBesselTransform(hr, r, k)
    #Sk = 1 + 2.0*rho/np.pi**3  * hk
    Sk =  1 + rho * hk
    #Theory
    Fig1.plot(r,hr+1, color=lc, marker='None', linewidth=1.0, label='eta=' + str(np.round(eta,2)))
    Fig2.plot(r,cr, color=lc, linestyle='-', marker='None', linewidth=3.0)
    Fig3.plot(k,Sk,color=lc, linewidth=3.0)
    
    #Numerical calculation
    if irho == 0:
        Flag, hr, cr, er, hk, Sk = OZ_Functions.OZSolver_MSA_Iterative(r,k,Ur,rho,sigma=1.0, kT=1.0, maxiter=10000, w_old_start = 0.50, tol=1e-10)
    else:
        Flag, hr, cr, er, hk, Sk = OZ_Functions.OZSolver_MSA_Iterative(r,k,Ur,rho,sigma=1.0, kT=1.0, maxiter=10000, w_old_start = 0.50, tol=1e-10, hr_guess=hr_guess, cr_guess=cr_guess)
    if Flag == 0:    
        Fig1.plot(r,hr+1, color=lc, marker='d', linestyle='None')
        Fig2.plot(r,cr, color=lc, linestyle='None', marker='d') #, linewidth=3.0)
        Fig3.plot(k,1 + rho * hk, color=lc,  marker='d', linestyle='None') #linewidth=3.0)
        cr_guess = cr #for the next rho
        hr_guess = hr #for the next rho
    else:
        break

Fig1.legend(loc='best')
Fig1.set_ylim(bottom=0.50)
#Fig2.legend(loc='best') 
#Fig3.legend(loc='best')

Fig.tight_layout()
plt.show()