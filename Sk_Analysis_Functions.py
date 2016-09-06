import numpy as np
from scipy.optimize import leastsq
import sys

def Error_SkPowerLaw(params, rho, Sk0, branch):
    Sk0_Pred = Fit_Sk0_PowerLaw(params, rho, branch)
    ErrorSq  = (Sk0_Pred/Sk0 - 1.0)**2 #relative error. 
    return ErrorSq

def Fit_Sk0_PowerLaw(params, rho, branch):
    c, rhostar, exponent = params
    if branch == 'right': #Liquid side
        Sk0_Pred = c / (rho - rhostar)**exponent
    elif branch == 'left':
        Sk0_Pred = c / (rhostar - rho)**exponent
    else:
        print "Branch should be left or right"
        sys.exit()
            

    return Sk0_Pred

def Fit(params_init_guess, rho, Sk0, branch):
    params_opt, pcov = leastsq(Error_SkPowerLaw, params_init_guess, args=(rho, Sk0, branch))
    return params_opt

