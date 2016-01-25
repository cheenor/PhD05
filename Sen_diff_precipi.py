#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 29 10:30:39 2015

@author: chenjh
"""
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
from pylab import *
import string
import numpy as np
import datetime
import scipy.stats as scista
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['contour.negative_linestyle'] = 'dashed'
mpl.rcParams['ytick.labelsize'] = 16
mpl.rcParams['xtick.labelsize'] = 16
plt.rc('lines', linewidth=4)
def readAscii(fpath,iskp,*nl):
    #iskp  the total line skipped of the file
    # fpath   the full path of the file
    # usage: onedim=readAscii(fpaht,iskp)
    onedim=[]
    linesplit=[]
    f=open(fpath)
    if nl:
        nrl=nl[0]
        ff=f.readlines()[iskp:nrl]  ## first line in obs file is legend 
    else:
        ff=f.readlines()[iskp:]
    for line in ff:
        line=string.lstrip(line)
        linesplit.append(line[:-1].split(' '))
    for lnstrs in linesplit:
        for strs in lnstrs:
#            if strs!='':
#                onedim.append(string.atof(strs))
            try:
                onedim.append(string.atof(strs))
            except:
                strs=""
#            else:
#                onedim.append(string.atof(strs))
    del linesplit,ff
    f.close()
    print len(onedim)
    return onedim
###############################################################################
MY3k=202
MY500=1202
nt=2880
CASE=['MLYRCTR_H','WTPCTR_H','ETPCTR_H']
Nsen=4
astr=[r'$(a)$',r'$(b)$', r'$(c)$',r'$(d)$',r'$(e)$',r'$(f)$'] 
DATESTR  =[[2010,6,2] , [2010,7,3] , [2010,6,3] ]
dirin500m='Z:/CRM/500m/'
dirinsen='Z:/CRM/SHLH/'
dirpic="D:/MyPaper/PhD05/Pics/"
dirobs='D:/MyPaper/PhD04/Data/TRMM/3B40/'
ydat_r=[ -50.000 ,    50.000 ,   164.286,    307.143,    478.571  ,  678.571 ,
      907.143 ,  1164.286,   1450.000,   1764.286 ,  2107.143,   2478.572 ,
      2878.572,   3307.143,  3764.286,  4250.000,   4764.286,   5307.143, 
      5878.571,   6478.571,   7107.143,  7764.286,  8450.000,  9164.285,  
      9907.143,  10678.570,  11478.570,  12307.143,  13164.285,  14050.000,
      14964.285,  15907.143,  16878.572,  17878.572,  18907.145,  19964.285,
      21050.000,  22164.285,  23307.145,  24478.572,  25678.572,  26907.145,
      28164.285,  29450.000,  30764.285,  32107.145,  33478.570,  34878.570,
      36307.141,  37764.285,  39250.000,  40750.000]
ydat=[]
for yd in ydat_r:
    ydat.append(yd*0.001)
nga=len(CASE)
pre3k=np.zeros(shape=(nga,MY3k,nt),dtype=float)
pbl3k=np.zeros(shape=(nga,MY3k,nt),dtype=float)
pre500=np.zeros(shape=(nga,MY500,nt),dtype=float)
pbl500=np.zeros(shape=(nga,MY500,nt),dtype=float)
presen=np.zeros(shape=(Nsen,nga,MY500,nt),dtype=float)
pblsen=np.zeros(shape=(Nsen,nga,MY500,nt),dtype=float)
alldate=[]
alltitle=[]
rainobs=[]
allregion=[]
for ig in range(0,nga):
    name500m=CASE[ig]
    cdate=DATESTR[ig]
    cdatestr='%4i'%cdate[0]+'%2.2i'%cdate[1]+'%2.2i'%cdate[2]
    if name500m[0:3]=='MLY' :
        region=name500m[0:4]     
    else:
        region=name500m[0:3]
    allregion.append(region)   
    fpath3k=dirin3k+region+'/'+fold3k+'/Simulation/'+'preci_'+name3k+'.txt'
    iskp=0
    onedim=readAscii(fpath3k,iskp)
    iskp=MY3k*4
    for it in range(0,nt):
        for ix in range(0,MY3k):
            ikk=it*iskp+MY3k*0+ix
            pre3k[ig,ix,it]=onedim[ikk]*1000.*3600
            ikk=it*iskp+MY3k*1+ix
            pbl3k[ig,ix,it]=onedim[ikk] 
    del onedim    
    fpath500m=dirin500m+region+'/run/postdata/'+'preci_'+name500m+'.txt'
    iskp=0
    onedim=readAscii(fpath500m,iskp)
    iskp=MY500*4
    for it in range(0,nt):
        for ix in range(0,MY500):
            ikk=it*iskp+MY500*0+ix
            pre500[ig,ix,it]=onedim[ikk]*1000.*3600 
            ikk=it*iskp+MY500*1+ix
            pbl500[ig,ix,it]=onedim[ikk] 
    !
    for isn in range(0,Nsen):
        fpathsen=dirinsen+region+'/run'+'%2.2i'%(isn+1)+'/postdata/'+'preci_'+name500m+'.txt'
        iskp=0
        onedim=readAscii(fpathsen,iskp)
        iskp=MY500*4
        for it in range(0,nt):
            for ix in range(0,MY500):
                ikk=it*iskp+MY500*0+ix
                presen[isn,ig,ix,it]=onedim[ikk]*1000.*3600 
                ikk=it*iskp+MY500*1+ix
                pblsen[isn,ig,ix,it]=onedim[ikk] 
    alltitle.append(astr[ig]+' '+region)
    year=cdate[0]
    month=cdate[1]
    day=cdate[2]
    datestart=datetime.datetime(year,month,day,0,0,0)
    det=datetime.timedelta(hours=3)            
    dateiso=[]            
    for dt in range(0,nt/12):
        dateiso.append(datestart+dt*det)
    xdate=[]              
    for tm in dateiso:
        xdate.append(datetime.datetime.strftime(tm,"%d/%b")) 
    alldate.append(xdate)
    fpath=dirobs+region+'-'+ cdatestr+'_031d_TRMM3B40.txt'
    iskp=1
    rain=readAscii(fpath,iskp)
    for i in range(0,nt/12):
        if rain[i]<0:
            rain[i]=0
    rainobs.append(rain)
###############################################################################
presenma=np.zeros(shape=(Nsen+1,nga,nt),dtype=float)
pblsenma=np.zeros(shape=(Nsen+1,nga,nt),dtype=float)
#pre500ma=np.zeros(shape=(nga,nt),dtype=float)
#pbl500ma=np.zeros(shape=(nga,nt),dtype=float)
presenma3h=np.zeros(shape=(Nsen+1,nga,nt/12),dtype=float)
pblsenma3h=np.zeros(shape=(Nsen+1,nga,nt/12),dtype=float)
#pre500ma3h=np.zeros(shape=(nga,nt/12),dtype=float)
#pbl500ma3h=np.zeros(shape=(nga,nt/12),dtype=float)
for ig in range(0,nga):
    for it in range(0,nt):
        for isn in range(1,Nsen+1):
            presenma[isn,ig,it]=presen[isn,ig,1:MY3k-1,it].mean()
            pblsenma[isn,ig,it]=pblsen[isn,ig,1:MY3k-1,it].mean()
        presenma[0,ig,it]=pre500[ig,1:MY500-1,it].mean()
        presenma[0,ig,it]=pbl500[ig,1:MY500-1,it].mean()
    for itt in range(0,nt/12):
        nt1=itt*12-6
        nt2=itt*12+6
        if itt==0:
            nt1=0
            nt2=7
        if itt==nt/12-1:
            nt1=nt-6
            nt2=nt
        dn=float(nt2-nt1)
        for n in range(nt1,nt2):
            for isn in range(0,Nsen+1):
                presenma3h[isn,ig,itt]= presenma3h[isn,ig,itt]+ presenma[ig,n]/dn
                pblsenma3h[isn,ig,itt]= pblsenma3h[isn,ig,itt]+ pblsenma[ig,n]/dn
                # if isn =0, it's 500m normal
#            pre500ma3h[ig,itt]= pre500ma3h[ig,itt]+ pre500ma[ig,n]/dn
#            pbl500ma3h[ig,itt]= pbl500ma3h[ig,itt]+ pbl500ma[ig,n]/dn
###############################################################################
xxx=range(0,ns)
for ig in range(0,nga):
    casename=allregion[iga]
    xdate=allxdate[iga]
    nr=int((Nsen+1.)/2)
    fig,ax=plt.subplots(nrows=2,ncols=3,figsize=(18,9))
    color_cycle=['deeppink', 'lime', 'b', 'y','indigo', 'cyan']
    wd=[2,2,2,2,2]
    jr=0
    jc=0
    ns=nt/12
    for isn in range(0,Nsen+1):
        if jc==3:
            jc=0
            jr=jr+1
        print jr,jc
        if labstr="Run"+"%2.2i"%isn
        ax[jr,jc].plot(xxx[0:ns],pre3kma3h[ig,0:ns],color=color_cycle[0],lw=wd[0],label=labstr) 
        ax[jr,jc].plot(xxx[0:ns],rainobs[ig][0:ns],color=color_cycle[1],lw=wd[1],label="TRMM")
        #ax[jr,jc].plot(xxx[0:ns],pre500ma3h[ig,0:ns],color=color_cycle[2],lw=wd[2],label="CRM 500m")
        ax[jr,jc].set_title(alltitle[ig], fontsize=18) 
        """      
        xx1=[]
        yy1=[]
        yy2=[]
        for irn in range(0,ns):
            if pre3kma3h[ig,irn]>0.0:
                xx1.append(pre3kma3h[ig,irn])
                yy1.append(rainobs[ig][irn])
                yy2.append(pre500ma3h[ig,irn])
        r1,p1=scista.pearsonr(xx1,yy1)
        r1str="r1= "+"%.2f"%r1
        r2,p2=scista.pearsonr(xx1,yy2)
        r2str="r2= "+"%.2f"%r2
        """
        if isn==5 :
            ax[jr,jc].legend(loc=(0.65,0.65),frameon=False)
        """
        if isn>3 :
            ax[jr,jc].set_yticks(range(0,3))
            ax[jr,jc].text(4,1.8,r1str,fontsize=18)
            ax[jr,jc].text(4,1.6,r2str,fontsize=18)
        else:
            ax[jr,jc].set_yticks(range(0,7))
            ax[jr,jc].text(4,5.4,r1str,fontsize=18)
            ax[jr,jc].text(4,4.8,r2str,fontsize=18)
        """
        ax[jr,jc].set_xticks(range(0,ns,72))
        xticklabels = [xdate[nn] for nn in range(0,ns,72)] 
        ax[jr,jc].set_xticklabels(xticklabels, size=18)
        jc=jc+1
    fig.text(0.03, 0.7, r'Precipitation ($mm$ $hr^{-1}$)', ha = 'left',fontsize=24,rotation=90)
    fig.subplots_adjust(left=0.1,bottom=0.12,right=1-0.05,top=1-0.1,hspace=0.4)
    plt.show()                     
    plt.savefig(dirpic+'Sen_'+casename+'_rainfall_Color.png',dpi=300)          
    plt.show()
    plt.close()