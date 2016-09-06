import OZ.Potentials as Potentials
import OZ.OZ_Functions as OZF
import OZ.PP_Functions  as PP_Functions
#import FigFuncs
#import HNC
import OZ.RHNC as RHNC


def SinglePointCalculation(r, k, Ur, Ur_ref, T, rho, OutFile=None, cr_guess=None):    

    #RHNC Solver
    Flag, hr, cr, er, hk, Sk = RHNC.OZSolver_RHNC_Iterative(r, k, Ur, Ur_ref, rho, kT=T, maxiter=10000, w_old_start=0.50,w_old_max=0.99,tol=1e-8, cr_guess=cr_guess)
    B2 = Potentials.CalcB2(r,Ur,kT=T)
    mu_ex, mu = PP_Functions.ExcessChemPot(r, hr, cr, rho,kT=T)

    rmin = 0.50
    if Flag == 0:

        gr = hr + 1.0
        
        P_virial = PP_Functions.Calc_Pvirial(r, gr, Ur, rho, kT=T, rmin=0.50)
        s2 = PP_Functions.s2(r, gr, rho)

        ListHeader=['T', 'rho', 'B2', 'P_virial' ,'Sk0', 'kappa', 'mu', 's2']
        ListValues = [T, rho, B2, P_virial, Sk[0], Sk[0]/(rho*T), mu, s2]

        

        if OutFile != None:
            OutFileH = open(OutFile, 'w')
            OutFileH.write(OZF.ListToTabbedStr(ListHeader))
            OutFileH.write(OZF.ListToTabbedStr(ListValues))

        return ListHeader, ListValues, gr, cr, Sk

    else:
        print "Not converged"
        return 0
