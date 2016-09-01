# -*- coding: utf-8 -*-
import OZ_Functions as OZF
import numpy as np 

def OZSolver_MSA(r,k,Ur,rho,sigma=1.0, kT=1.0, maxiter=10000, w_old = 0.80, tol=1e-8, hr_guess=None, cr_guess=None):
    '''
    Ornstein-Zerinike solver with MSA closure for hard core + long range attraction/repulsion. 
    '''
    nr = len(r)
    w_new = 1.0 - w_old
    #Initial guess for hr/hk    
    if hr_guess == None:
        hr = np.where(r < sigma, -1, np.exp(-Ur/kT)-1) #Ideal gas (guess) + hard core inside sigma.
    else:
        hr = hr_guess    
    hk = OZF.FourierBesselTransform(hr,r,k)
    #Initial guess for c(r)/ck: c(r)  = 0 outside hard core and inside, we initially guess as 0. 
    if cr_guess == None:
        cr = np.where(r > sigma, -Ur/kT, 0)
    else:
        cr = cr_guess
    ck = OZF.FourierBesselTransform(cr, r, k)
    
    residual = 1e10
    iter = 0
    #hk_new = hk
    while residual > tol and iter < maxiter:        
        ck_temp = hk/(1.0 + rho * hk) #Ornstein-Zernike relation
        cr_temp = OZF.InvFourierBesselTransform(ck_temp,r,k)
        cr_new = np.where(r > sigma, -Ur/kT, cr_temp) #enforcing c(r) = exp(-beta * U(r)) outside the hard core. 
                
        hk_temp = ck/(1.0 - rho * ck)    
        hr_temp = OZF.InvFourierBesselTransform(hk_temp, r, k)
        hr_new = np.where(r < sigma, -1, hr_temp) #enforcing gr(r) = 0 within the hard core.         
                
        err = np.where(r > sigma, hr_new - hr, cr_new-cr)
        residual_new = np.linalg.norm(err)
        if residual_new < tol:
            break 
        else:
            hr = hr_new * w_new + hr * w_old
            cr = cr_new * w_new + cr * w_old
            
            hk = OZF.FourierBesselTransform(hr,r,k)
            ck = OZF.FourierBesselTransform(cr,r,k)
           
            residual = residual_new
            iter = iter + 1
            #print "iter = ", iter, "residual = ", residual
        
    if residual_new != np.inf and residual_new < tol:
        hk = OZF.FourierBesselTransform(hr,r,k)        
        Sk = 1 + rho * hk
        #check for long range oscillatory behavior in g(r) or negative values of S(k)
        if np.mean(abs(hr[int(nr/2):])) > 1e-4 or np.min(Sk) < 0.0:
            Flag = -1
            print "OZSolver: Converged, but UNPHYSICAL result - S(k) < 0.0"
        else:
            Flag = 0
            print "OZSolver: Converged, and physical result - S(k) > 0.0"
            print "rho = ", rho, "w_old = ", w_old
                        
    else:
        print "Not converged to below tolerance"
        print "Residual = ", residual_new
        Flag = -1

    if Flag  == 0:
        return Flag, hr, cr, hr-cr, hk, Sk
    else:        
        return Flag , np.zeros(nr), np.zeros(nr), np.zeros(nr), np.zeros(nr), np.zeros(nr)


def OZSolver_MSA_Iterative(r, k, Ur, rho, sigma=1.0, kT=1.0, maxiter_OZ = 20, maxiter=10000, w_old_start=0.50, w_old_max=0.99,tol=1e-8,hr_guess=None, cr_guess=None):
    '''
    Solve the OZ equation by changing the increasing the weight to
    picard iterator in case of non-convergence for small values of w_old
    '''
    Flag = -1
    w = w_old_start
    iter_OZ = 0
    while Flag == -1 and w < w_old_max * 0.95:
        iter_OZ = iter_OZ + 1
        Flag, hr, cr, er, hk, Sk = OZSolver_MSA(r,k,Ur,rho,sigma=1.0, kT=1.0, maxiter=10000, w_old = 0.80, tol=1e-8, hr_guess=hr_guess, cr_guess=cr_guess)
        print "iter_OZ = ", iter_OZ, "Setting w_old = ", w, "Flag = ", Flag, "rho = ", rho
            
        if Flag == 0:
            break            
        else:
            w = w + (w_old_max-w) * 0.10  #update picard weight.
            
    if Flag == 0:
        return Flag, hr, cr, er, hk, Sk
    else:
        return Flag, np.zeros_like(r), np.zeros_like(r), np.zeros_like(r), np.zeros_like(r), np.zeros_like(r)
