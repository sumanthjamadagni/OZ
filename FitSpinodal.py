import numpy as np
import OZ.Sk_Analysis_Functions as SkF
import sys
import matplotlib.pyplot as plt
import matplotlib
from itertools import cycle

colors= ['red', 'green', 'black', 'cyan', 'blue'] # 'magenta']
colorcycler = cycle(colors)

 

FS = 18
matplotlib.rc('xtick', labelsize=FS)
matplotlib.rc('ytick', labelsize=FS)

T_Str = ['0.80', '0.90', '1.00', '1.00', '1.10', '1.20', '1.30']
branch='right'
Spinodal = []
T_float = []

Fig = plt.figure(figsize=[15,6])
Fig1 = Fig.add_subplot(121)
Fig1.set_ylabel('Sk[0]', fontsize=FS)
Fig1.set_xlabel('$\\rho$', fontsize=FS)


Fig2 = Fig.add_subplot(122)
Fig2.set_xlabel('$\\rho_{spinodal}$', fontsize=FS)
Fig2.set_ylabel('T', fontsize=FS)


for T in T_Str:
    c = next(colorcycler)
    print "c = ", c
    FileName = 'T' + T + '-Liq.dat'
    Data = np.loadtxt(FileName, skiprows=1)
    rho = Data[:,1]
    Sk0 = Data[:,4]
    #initial guess for rho_star slightly smaller than the smallest density for which the calculation converged
    params_init_guess = (1.0, np.min(rho)-0.05, 1.5)            
    params_opt =  SkF.Fit(params_init_guess, rho, Sk0, branch)

    c, rho_sp, exponent = params_opt
    Spinodal.append(rho_sp)
    T_float.append(float(T))
    
    Sk0_Fit = SkF.Fit_Sk0_PowerLaw(params_opt, rho, branch)

    Fig1.plot(rho, Sk0, linestyle='None', linewidth=2, marker='d', label = T + ',liq')
    Fig1.plot(rho, Sk0, linestyle='-', linewidth=2)

    print T, params_opt

Fig2.plot(np.array(Spinodal), np.array(T_float), 'r-d')

#-------------------- REPEAT FOR VAPOR SIDE
Spinodal = []
T_float = []
branch='left'
for T in T_Str:
    c = next(colorcycler)
    print "c = ", c
    FileName = 'T' + T + '-Vap.dat'
    Data = np.loadtxt(FileName, skiprows=1)
    rho = Data[:,1]
    Sk0 = Data[:,4]
    #initial guess for rho_star slightly larger than the largest density for which the calculation converged
    params_init_guess = (1.0, np.max(rho)*1.25, 1.5)            
    params_opt =  SkF.Fit(params_init_guess, rho, Sk0, branch)

    c, rho_sp, exponent = params_opt
    Spinodal.append(rho_sp)
    T_float.append(float(T))
    
    Sk0_Fit = SkF.Fit_Sk0_PowerLaw(params_opt, rho, branch)

    Fig1.plot(rho, Sk0, linestyle='None', linewidth=2, marker='d', label = T + ',vap')
    Fig1.plot(rho, Sk0_Fit, linestyle='-', linewidth=2)

    print T, params_opt

Fig2.plot(np.array(Spinodal), np.array(T_float), 'b-d')

Fig1.legend(loc='best', ncol=2)





Fig.tight_layout()
plt.show()

Fig.savefig('LJ_Spinodal.pdf')
