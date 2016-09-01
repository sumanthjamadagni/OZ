# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib

#colors = ['red', 'blue', 'green', 'cyan', 'black', 'orange']
#ColorCycler = cycle(colors)

FS = 18
matplotlib.rc('xtick', labelsize=FS)
matplotlib.rc('ytick', labelsize=FS)


def CreateFig(x,y,title=None, xlabel=None, ylabel=None, xleft=None, xright=None, ybottom=None, ytop=None, lw=2, marker='d', ls='-'):
    Fig = plt.figure()
    Fig1 = Fig.add_subplot(111)
    print "xright = ", xright
    if xlabel:
        Fig1.set_xlabel(xlabel, fontsize=FS)

    if ylabel:
        Fig1.set_ylabel(ylabel, fontsize=FS)

    if title:
        Fig1.set_title(title, fontsize=FS)

    Fig1.plot(x, y, linewidth=lw, linestyle=ls, marker=marker)
    if xleft is not None and xright is not None:    
        Fig1.set_xlim([xleft,xright])
    if ybottom is not None and ytop is not None:    
        Fig1.set_ylim([ybottom, ytop])
        
    return Fig
    
def CreateFigObject(title=None, xlabel=None, ylabel=None):
    Fig = plt.figure()
    Fig1 = Fig.add_subplot(111)
    if xlabel:
        Fig1.set_xlabel(xlabel)

    if ylabel:
        Fig1.set_ylabel(ylabel)

    if title:
        Fig1.set_title(title)
        
    return Fig
    
