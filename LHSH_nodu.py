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
picout='D:/MyPaper/PhD05/Pics/'
Regions=["PRD","MLYR","NPC","NEC","WTP","ETP"]
datestrs=['20120401' , '20100602' , '20100802' ,
         '20120706' , '20100703' , '20100603' ]
n=len(Regions)
nn=249
for i in range(0,n):
    area=Regions[i]
    datestr=datestrs[i]
    iy=string.atoi(datestr[0:4])
    im=string.atoi(datestr[4:6])
    jd=string.atoi(datestr[6:8])
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
    #
    alldatare2=np.zeros(shape=(4,nn),dtype=float)
    for iv in range(0,4):
        for it in range(0,nn):
            its=it-4
            ite=it+4
            if its<0:
                its=0
                ite=ite-its
            if ite>nn:
                its=its-(nn-ite)
                ite=nn
            tmp=0.
            for j in range(its,ite):
                tmp=tmp+alldata[iv,j]
            alldatare2[iv,it]=tmp/9.
    fpath=dirin+area+'/'+area+'_'+datestr+'_031d_SHLH_ERA_nodu2.43'
    fout=open(fpath,'w')    
    for it in range(0,nn):
        for iv in range(0,4):
            if iv>1:
               alldatare2[iv,it]=alldata[iv,it] 
            outstr='%e '%alldatare2[iv,it]
            fout.write(outstr)
        fout.write('\n')     
    #######################################################################
    datestart=datetime.datetime(iy,im,jd,0,0,0)
    det=datetime.timedelta(hours=3)            
    dateiso=[] 
    xdate=[]                
    for dt in range(0,nn):
        dateiso.append(datestart+dt*det)           
    for tm in dateiso:
        xdate.append(datetime.datetime.strftime(tm,"%d/%b"))
    xxx=range(0,nn)
    charsize=16       
    fig,(ax0,ax1,ax2) = plt.subplots(nrows=3,ncols=1,figsize=(18,6))
    ax0.plot(xxx,alldata[0,:],'g',label='0')#[0:lcc-3])
    #plt.axis([0, nn, -100, 500])  ## x axis  y axis
    text1=r"($a$) Normal"
    ax0.text(2,203,text1,fontsize=16)
    ax0.set_xticks(range(0,nn,16))
    xticklabels = [xdate[i] for i in range(0,nn,16)] 
    ax0.set_xticklabels(xticklabels, size=charsize)
    #
    ax1.plot(xxx,alldatare[0,:],'g',label='1')#[0:lcc-3])
    #plt.axis([0, nn, -100, 500])  ## x axis  y axis
    text1=r"($b$) Day Mean"
    ax1.text(2,203,text1,fontsize=16)
    ax1.set_xticks(range(0,nn,16))
    xticklabels = [xdate[i] for i in range(0,nn,16)] 
    ax1.set_xticklabels(xticklabels, size=charsize)
    #
    ax2.plot(xxx,alldatare2[0,:],'g',label='0')#[0:lcc-3])
    #plt.axis([0, nn, -100, 500])  ## x axis  y axis
    text1=r"($c$) Moving average"
    ax2.text(2,203,text1,fontsize=16)
    ax2.set_xticks(range(0,nn,16))
    xticklabels = [xdate[i] for i in range(0,nn,16)] 
    ax2.set_xticklabels(xticklabels, size=charsize)
    plt.show()                     
    plt.savefig(picout+area+'_surface.png',dpi=300)          
    plt.show()
    plt.close()        
        
            
