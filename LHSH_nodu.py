# -*- coding: utf-8 -*-
"""
Created on Sat Dec 05 14:50:20 2015

@author: chenjh
"""
import matplotlib as mpl
import numpy as np
import matplotlib.cm as cm
import datetime
from pylab import *
import matplotlib.pyplot as plt
import matplotlib.dates as matdate
import calendar
import string
import time
matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'
matplotlib.rcParams['contour.negative_linestyle'] = 'dashed'
matplotlib.rcParams['ytick.labelsize'] = 16
matplotlib.rcParams['xtick.labelsize'] = 16
def readAscii(fpath,iskp):
    #iskp  the total line skipped of the file
    # fpath   the full path of the file
    # usage: onedim=readAscii(fpaht,iskp)
    onedim=[]
    linesplit=[]
    f=open(fpath)
    ff=f.readlines()[iskp:]  ## first line in obs file is legend 
    for line in ff:
        line=string.lstrip(line)
        linesplit.append(line[:-1].split(' '))
    for lnstrs in linesplit:
        for strs in lnstrs:
            if strs!='':
                onedim.append(string.atof(strs))
    del linesplit
    f.close()
    return onedim
###############################################################################
dirin='D:/MyPaper/PhD04/Cases/ERA/FORCING/'
Regions=["PRD","MLYR","NPC","NEC","WTP","ETP"]
datestrs=['20120401' , '20100602' , '20100802' ,
         '20120706' , '20100703' , '20100603' ]
n=len(Regions)
nn=249
for i in range(0,n):
    area=Regions[i]
    datestr=datestrs[i]
    fpath=dirin+area+'/'+area+'_'+datestr+'_031d_SHLH_ERA.43'
    iskp=0
    onedim=readAscii(fpath,iskp)
    alldata=np.zeros(shape=(4,nn),dtype=float)
    for it in range(0,nn):
        for iv in range(0,4):
            k=it*4+iv
            alldata[iv,it]=onedim[k]
    alldatare=np.zeros(shape=(4,nn),dtype=float)    
    for iv in range(0,4):
        tmp=0.
        for j in range(0,2):
            tmp=tmp+alldata[iv,j]
        for j in range(0,2):        
            alldatare[iv,j]=tmp/2.
        for it in range(2,nn-7,8):
            tmp=0.            
            for j in range(0,8):
                k=it+j
                tmp=tmp+alldata[iv,k]
            for j in range(0,8):
                kk=it+j
                alldatare[iv,kk]=tmp/8.
        tmp=0.
        for j in range(nn-7,nn):
            tmp=tmp+alldata[iv,j]
        for j in range(nn-7,nn):        
            alldatare[iv,j]=tmp/7.
    fpath=dirin+area+'/'+area+'_'+datestr+'_031d_SHLH_ERA_nodu.43'
    fout=open(fpath,'w')    
    for it in range(0,nn):
        for iv in range(0,4):
            if iv>1:
               alldatare[iv,it]=alldata[iv,it] 
            outstr='%e '%alldatare[iv,it]
            fout.write(outstr)
        fout.write('\n')     
            
            
        
            
