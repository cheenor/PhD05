# -*- coding: utf-8 -*-
"""
Created on Fri Dec 25 15:51:27 2015

@author: chenjh
"""
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
from pylab import *
import string
import numpy as np
import datetime
import struct
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['contour.negative_linestyle'] = 'dashed'
mpl.rcParams['ytick.labelsize'] = 16
mpl.rcParams['xtick.labelsize'] = 16
plt.rc('lines', linewidth=4)
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
def readBin(fpath):
    # fpath   the full path of the file
    # usage: onedim=readAscii(fpath)
    onedim=[]
    binfile=open(fpath,'rb')
    bins=binfile.read()
    nn= len(bins)/4
    for n in range(0,nn):
        k1=n*4
        k2=k1+4
        aa,=struct.unpack('f',bins[k1:k2])
        onedim.append(aa)        
    return onedim
def getprof(x):
    nz=len(x[:,0])
    nx=len(x[0,:])
    bk=np.zeros(shape=(nz),dtype=float)
    for iz in range(0,nz):
        bk[iz]=0.
        for n in range(0,nx):
            bk[iz]=bk[iz]+x[iz,n]/(nx*1.0)
    return bk
###############################################################################
nt=121
nday=31
nz=52
CASE=['PRDCTR_H','MLYRCTR_H','NPCCTR_H','NECCTR_H','WTPCTR_H','ETPCTR_H']
astr=[r'$(a)$',r'$(b)$', r'$(c)$',r'$(d)$',r'$(e)$',r'$(f)$'] 
DATESTR  =[[2012,4,1] , [2010,6,2] , [2010,8,2] ,
           [2012,7,6] , [2010,7,3] , [2010,6,3] ]
nx=1202
nga=len(CASE)
dirin3k='D:/MyPaper/PhD04/Cases/'
dirin500m='Z:/CRM/500m/'
dirpic='D:/MyPaper/PhD05/Pics/'
#
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
xdat=range(0,nt)  
###############################################################################
"""    
smt=np.ndarray(shape=(nz,nt,nga),dtype=float)
smq=np.ndarray(shape=(nz,nt,nga),dtype=float)
smu=np.ndarray(shape=(nz,nt,nga),dtype=float)
smw=np.ndarray(shape=(nz,nt,nga),dtype=float)
smqc=np.ndarray(shape=(nz,nt,nga),dtype=float)
smqr=np.ndarray(shape=(nz,nt,nga),dtype=float)
smqa=np.ndarray(shape=(nz,nt,nga),dtype=float)
smqb=np.ndarray(shape=(nz,nt,nga),dtype=float)
smrh=np.ndarray(shape=(nz,nt,nga),dtype=float)
smfqv=np.ndarray(shape=(nz,nt,nga),dtype=float)
smftc=np.ndarray(shape=(nz,nt,nga),dtype=float)
smrat=np.ndarray(shape=(nz,nt,nga),dtype=float)
smtc=np.ndarray(shape=(nz,nt,nga),dtype=float)
"""
sm3k=np.ndarray(shape=(13,nz,nt,nga),dtype=float)
sm500m=np.ndarray(shape=(13,nz,nt,nga),dtype=float)
alltitle=[]
allxdate=[]
for iga in range(0,nga):
    name500m=CASE[iga]
    cdate=DATESTR[iga]
    cdatestr='%4i'%cdate[0]+'%2.2i'%cdate[1]+'%2.2i'%cdate[2]
    if name500m[0:3]=='MLY' :
        region=name500m[0:4]
        name3k=name500m[0:7]+'_EC'       
    else:
        region=name500m[0:3]
        name3k=name500m[0:6]+'_EC'
    fold3k='CTREC'+cdatestr    
    fpath3k=dirin3k+region+'/'+fold3k+'/Simulation/'+name3k+'_All.txt'
    iskp=0
    onedim=readAscii(fpath3k,iskp)
    iskp=13*nz
    ikkk=0
    for it in range(0,nt):
        for iz in range(0,nz):
            for ik in range(0,13):
                ikk=it*iskp+nz*ik+iz
                sm3k[ik,iz,it,iga]=onedim[ikk]   
    del onedim    
    fpath500m=dirin500m+region+'/run/'+name500m+'_All.txt'
    iskp=0
    onedim=readAscii(fpath500m,iskp)
    iskp=13*nz
    ikkk=0
    for it in range(0,nt):
        for iz in range(0,nz):
            for ik in range(0,13):
                ikk=it*iskp+nz*ik+iz
                sm500m[ik,iz,it,iga]=onedim[ikk] 
    alltitle.append(astr[iga]+' '+region)
    year=cdate[0]
    month=cdate[1]
    day=cdate[2]
    datestart=datetime.datetime(year,month,day,0,0,0)
    det=datetime.timedelta(hours=6)            
    dateiso=[]            
    for dt in range(0,nt):
        dateiso.append(datestart+dt*det)
    xdate=[]              
    for tm in dateiso:
        xdate.append(datetime.datetime.strftime(tm,"%d/%b")) 
    allxdate.append(xdate)
###############################################################################
fig,ax=plt.subplots(nrows=3,ncols=2,figsize=(18,12))
#fig,axs=plt.subplots(nrows=2,ncols=3,figsize=(12,12))
color_cycle=['deeppink', 'lime', 'b', 'y','indigo', 'cyan']
wd=[2,2,2,2,2]
jc=0
jr=0
ij=1
for iga in range(0,nga):
    if jc==2:
        jc=0
        jr=jr+1
    mker=alltitle[iga]
    xdate=allxdate[iga]
    levs1=[-15,-12,-6,-3,4,8,12,15]
    colors1=['g','g','g','g','r','r','r','r']
    linetype1=['dotted','dotted','dotted','dotted','solid','solid','solid','solid'] 
    levs2=[-3,-1,1,2,3,4]
    colors2=['g','g','r','r','r','r']
    linetype2=['dotted','dotted','solid','solid','solid','solid']
    titlename=[r"Temperature Bias ($K$)",r"Water Vapor Mixing Ratio Bias ($g$ $kg^{-1}$) "]
    font = {'family' : 'serif',
        'color'  : 'k',
        'weight' : 'normal',
        'size'   : 16,
        }     
    plt.subplot(3,2,ij)
    zdat=sm500m[0,:,:,iga]-sm3k[0,:,:,iga]
    #zdat=smtco   #-obstmp
    ##zdat[0,:]=0.0   ## the first level is below surface ground
    ax[jr,jc]=plt.contour(xdat,ydat,zdat,colors='r',
        linewidths=1.5,levels=levs1,linestyles=linetype1)                           
#    plt.title(titlename[0],fontsize=16)                          
    plt.axis([0, 121, 0, 16])
    plt.clabel(ax[jr,jc],inline=1,fmt='%1d',fontsize=12)
    zdat=sm500m[1,:,:,iga]-sm3k[1,:,:,iga]
    zdat[0,:]=0.0   ## the first level is below surface ground
    ax[jr,jc]=plt.contour(xdat,ydat,zdat,colors='g',
        linewidths=1.5,levels=levs2,linestyles=linetype2)  
    plt.axis([0, 121, 0, 16])
    plt.clabel(ax[jr,jc],inline=1,fmt='%1d',fontsize=12)
    axx=fig.add_subplot(3,2,ij)
    text1=mker  #r"($a$)"
    axx.text(1.5,16.5,text1,fontsize=18)                        
    axx.set_xticks(range(0,nt,16))
    xticklabels = [xdate[nn] for nn in range(0,nt,16)] 
    axx.set_xticklabels(xticklabels, size=16)
    plt.ylabel(r'Height ($km$)', fontdict=font)
    jc=jc+1
    ij=ij+1                
plt.show()
fig.subplots_adjust(left=0.1,bottom=0.1,right=1-0.1,top=1-0.1,hspace=0.4)
plt.savefig(dirpic+"ALLCASE_T&qv_resulo_diff.png",dpi=300)          
plt.show()
plt.close()                
###############################################################################
fig,ax=plt.subplots(nrows=3,ncols=2,figsize=(18,12))
#fig,axs=plt.subplots(nrows=2,ncols=3,figsize=(12,12))
color_cycle=['deeppink', 'lime', 'b', 'y','indigo', 'cyan']
wd=[2,2,2,2,2]
jc=0
jr=0
ij=1
for iga in range(0,nga):
    if jc==2:
        jc=0
        jr=jr+1
    mker=alltitle[iga]
    xdate=allxdate[iga]
    levs1=[-15,-12,-6,-3,4,8,12,15]
    colors1=['g','g','g','g','r','r','r','r']
    linetype1=['dotted','dotted','dotted','dotted','solid','solid','solid','solid'] 
    levs2=[-3,-1,1,2,3,4]
    colors2=['g','g','r','r','r','r']
    linetype2=['dotted','dotted','solid','solid','solid','solid']
    titlename=[r"Temperature Bias ($K$)",r"Water Vapor Mixing Ratio Bias ($g$ $kg^{-1}$) "]
    font = {'family' : 'serif',
        'color'  : 'k',
        'weight' : 'normal',
        'size'   : 16,
        }     
    plt.subplot(3,2,ij)
    zdat1=sm500m[4,:,:,iga]-sm3k[4,:,:,iga]
    zdat2=sm500m[5,:,:,iga]-sm3k[5,:,:,iga]
    zdat3=sm500m[6,:,:,iga]-sm3k[6,:,:,iga]
    zdat4=sm500m[7,:,:,iga]-sm3k[7,:,:,iga]
    zdat=zdat1+zdat2
    #zdat=smtco   #-obstmp
    ##zdat[0,:]=0.0   ## the first level is below surface ground
    ax[jr,jc]=plt.contour(xdat,ydat,zdat,colors='r',
        linewidths=1.5,levels=levs1,linestyles=linetype1)                           
#    plt.title(titlename[0],fontsize=16)                          
    plt.axis([0, 121, 0, 16])
    plt.clabel(ax[jr,jc],inline=1,fmt='%1d',fontsize=12)
    zdat=zdat3+zdat4
    zdat[0,:]=0.0   ## the first level is below surface ground
    ax[jr,jc]=plt.contour(xdat,ydat,zdat,colors='g',
        linewidths=1.5,levels=levs2,linestyles=linetype2)  
    plt.axis([0, 121, 0, 16])
    plt.clabel(ax[jr,jc],inline=1,fmt='%1d',fontsize=12)
    axx=fig.add_subplot(3,2,ij)
    text1=mker  #r"($a$)"
    axx.text(1.5,16.5,text1,fontsize=18)                        
    axx.set_xticks(range(0,nt,16))
    xticklabels = [xdate[nn] for nn in range(0,nt,16)] 
    axx.set_xticklabels(xticklabels, size=16)
    plt.ylabel(r'Height ($km$)', fontdict=font)
    jc=jc+1
    ij=ij+1                
plt.show()
fig.subplots_adjust(left=0.1,bottom=0.1,right=1-0.1,top=1-0.1,hspace=0.4)
plt.savefig(dirpic+"ALLCASE_cloudwater_resulo_diff.png",dpi=300)          
plt.show()
plt.close()   
###############################################################################
fig,ax=plt.subplots(nrows=2,ncols=3,figsize=(10,10))
#fig,axs=plt.subplots(nrows=2,ncols=3,figsize=(12,12))
color_cycle=['deeppink', 'g', 'b', 'indigo','yellow', 'cyan']
lss=['-',':','--','-.','-']
wd=[3,3,3,3,3]
jc=0
jr=0
for iga in range(0,nga):
    if jc==3:
        jc=0
        jr=jr+1
    tmp1=sm500m[4,:,:,iga]-sm3k[4,:,:,iga]
    tmp2=sm500m[5,:,:,iga]-sm3k[5,:,:,iga]
    tmp3=sm500m[6,:,:,iga]-sm3k[6,:,:,iga]
    tmp4=sm500m[7,:,:,iga]-sm3k[7,:,:,iga]
    zdat1=getprof(tmp1)
    zdat2=getprof(tmp2)
    zdat3=getprof(tmp3)
    zdat4=getprof(tmp4)
    zdat5=zdat1+zdat2+zdat3+zdat4
    marknm=alltitle[iga]
    #zdat=smtco   #-obstmp
    ##zdat[0,:]=0.0   ## the first level is below surface ground
    ax[jr,jc].plot(zdat1,ydat,c=color_cycle[0], ls=lss[0],
        lw=3,label=r'qc')                                                    
    ax[jr,jc].plot(zdat2,ydat,c=color_cycle[1],ls=lss[1],
        lw=3,label=r'qr') 
    ax[jr,jc].plot(zdat3,ydat,c=color_cycle[2],ls=lss[2],
        lw=3,label=r'qa')         
    ax[jr,jc].plot(zdat4,ydat,c=color_cycle[3],ls=lss[3],
        lw=3,label=r'qb') 
    ax[jr,jc].plot(zdat5,ydat,c='k',ls=lss[4],
        lw=3,label=r'total') 
    if jc==0:    
        ax[jr,jc].set_ylabel(r'Height ($km$)', fontdict=font)
    if jc==2 and jr==1:
        ax[jr,jc].legend(frameon=False)
    ax[jr,jc].set_title(marknm,fontsize=14)
    ax[jr,jc].set_ylim(0,16)
    ax[jr,jc].set_xlim(-0.04,0.04)
    xmajorLocator   = MultipleLocator(0.015)
    ax[jr,jc].xaxis.set_major_locator(xmajorLocator)
    ymajorLocator   = MultipleLocator(4)
    ax[jr,jc].yaxis.set_major_locator(ymajorLocator)
    jc=jc+1        
plt.show()
fig.subplots_adjust(left=0.1,bottom=0.1,right=1-0.1,top=1-0.1,hspace=0.4)
plt.savefig(dirpic+"ALLCASE_cloudwater_resulo_diff_profile.png",dpi=300)          
plt.show()
plt.close()                 
                
                