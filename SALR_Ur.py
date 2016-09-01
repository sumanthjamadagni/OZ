# -*- coding: utf-8 -*-

# -*- coding: utf-8 -*-

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
D  = [2.0]
FigUr = plt.figure()
FigUr1 = FigUr.add_subplot(111)
#FigUr2 = FigUr.add_subplot(122)

FigUr1.axhline(0, linewidth=2, linestyle='--', color='black')
FigUr1.set_xlim([0.0,5])
d = 2
umin = 1e4
E = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
for e in E:
    c = next(ColorCycler)
    #UrR = Potentials.SRLRPotential(r,sig=1.0, eps=3.0, A=0.70, d=d, m=36, n=18)
    UrA = Potentials.SALRPotential(r,sig=1.0, eps=e, A=0.40, d=d, m=100, n=50)
    #FigUr1.plot(r,UrR,linewidth=2.0, label='d_R=' + str(d), color=c, linestyle='--')
    FigUr1.plot(r,UrA,linewidth=2.0, label='d_A=' + str(d), color=c, linestyle='-')
    umin = min(np.min(UrA), umin)

    T = 1.0    
    B2_Att = Potentials.CalcB2(r,UrA,kT=T)   
    print "e = ", e, "B2_Att = ", B2_Att #, "B2_Rep = ", B2_Rep


FigUr1.set_ylim([umin-1,4])
FigUr1.legend(loc='best')
plt.show()
