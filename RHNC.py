# -*- coding: utf-8 -*-
import numpy as np
import OZ_Functions as OZF
import HNC

def RHNC_Closure(er, Ur, Ur_ref, gr_ref, er_ref, kT=1.0):
    deltaU = Ur - Ur_ref
    delta_er = er - er_ref #difference in the indirect correlation function
    cr = gr_ref * np.exp(-deltaU/kT + delta_er) - er - 1.0        
    return cr

def OZSolver_RHNC(r, k, Ur, Ur_ref, rho, kT=1.0, maxiter=10000, w_old=0.50, tol=1e-10, cr_guess=None):     
    Flag_ref, hr_ref, cr_ref, er_ref, hk_ref, Sk_ref = HNC.OZSolver_HNC_Iterative(r, k, Ur_ref, rho, kT=1.0, maxiter=10000, w_old_start=0.50, tol=tol) 
    
    nr = len(r)
    if Flag_ref != 0:
        print "reference fluid HNC calculation failed!"
        return Flag_ref , np.zeros(nr), np.zeros(nr), np.zeros(nr), np.zeros(nr), np.zeros(nr)
    
                
    if Flag_ref == 0: #the HNC calculation for the refernce fluid has converged. 
        gr_ref = hr_ref + 1        

        w_new = 1.0 - w_old
        if cr_guess == None:
            cr = OZF.cr_initguess(Ur)
        else:
            print "Using initial guess for c(r) in HNC solver"
            cr = cr_guess 
            
        iter = 0        
        while iter < maxiter:    
            ck = OZF.FourierBesselTransform(cr,r,k)            
            ek = OZF.Calculate_ek(ck,rho)            
            er = OZF.InvFourierBesselTransform(ek, r, k)            
            cr_new = RHNC_Closure(er, Ur, Ur_ref, gr_ref, er_ref, kT=kT)

    
            residual = np.linalg.norm(cr-cr_new)
        
            if residual < tol:
                print "OZSolver: Converged!! Exiting loop"
                iter = 2 * maxiter  #to break the loop
            elif residual == np.inf:
                print "COZSolver: crashed!! Exiting loop"
                iter = 2 * maxiter  #to break the loop
            else: 
                #print "iteration  = ", iter, "residual = ", residual_new        
                cr = cr_new*w_new + cr*w_old #mixing old + new

            iter = iter + 1                
        
        hr = OZF.Calculate_hr(er, cr)
        #h(r) may have converged to non-physical results: 
        if residual != np.inf and residual < tol:
            hk = OZF.FourierBesselTransform(hr,r,k)            
            Sk = 1 + rho * hk
        #check for long range oscillatory behavior in g(r) or negative values of S(k)
            if np.mean(abs(hr[int(nr/2):])) > 1e-3 or np.min(Sk) < 0.0:
                Flag = -1
                print "OZSolver: Converged, but UNPHYSICAL result - S(k) < 0.0"
            else:
                Flag = 0
                print "OZSolver: Converged, and physical result - S(k) > 0.0"
                print "rho = ", rho, "w_old = ", w_old
            
            
        else:
            print "Not converged to below tolerance"
            print "Residual = ", residual
            Flag = -1

    if Flag  == 0:
        return Flag, hr, cr, er, hk, Sk
    else:        
        return Flag , np.zeros(nr), np.zeros(nr), np.zeros(nr), np.zeros(nr), np.zeros(nr)

def OZSolver_RHNC_Iterative(r, k, Ur, Ur_ref, rho, kT=1.0, maxiter=10000, w_old_start=0.50,w_old_max=0.99,tol=1e-8, cr_guess=None):
    '''
    Solve the OZ equation by changing the increasing the weight to
    picard iterator in case of non-convergence for small values of w_old
    '''

    Flag = -1
    w = w_old_start
    iter_OZ = 0
    while Flag == -1 and w < w_old_max * 0.95:
        iter_OZ = iter_OZ + 1
        
        Flag, hr, cr, er, hk, Sk = OZSolver_RHNC(r, k, Ur, Ur_ref, rho, kT=kT, maxiter=maxiter, w_old=w, tol=tol, cr_guess=cr_guess)
        print "iter_OZ = ", iter_OZ, "Setting w_old = ", w, "Flag = ", Flag, "rho = ", rho
            
        if Flag == 0:
            break            
        else:
            w = w + (w_old_max-w) * 0.10  #update picard weight.
            
            
    if Flag == 0:
        return Flag, hr, cr, er, hk, Sk
    else:
        return Flag, np.zeros_like(r), np.zeros_like(r), np.zeros_like(r), np.zeros_like(r), np.zeros_like(r)
    
