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
nt=2880
nday=31
nz=34
CASE=['MLYRCTR_H','WTPCTR_H','ETPCTR_H']
Nsen=4
astr=[r'$(a)$',r'$(b)$', r'$(c)$',r'$(d)$',r'$(e)$',r'$(f)$'] 
DATESTR  =[[2010,6,2] , [2010,7,3] , [2010,6,3] ]
dirin500m='Z:/CRM/500m/postdata/'
dirinsen='Z:/CRM/SHLH/postdata/'
dirpic="D:/MyPaper/PhD05/Pics/"
nga=len(CASE)
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
rawdatasen=np.zeros(shape=(Nsen+1,nga,2,nz,nt),dtype=float)
cloudpathsen=np.zeros(shape=(Nsen+1,nga,5,2,nt),dtype=float)
#rawdata500=np.zeros(shape=(nga,2,nz,nt),dtype=float)
#cloudpath500=np.zeros(shape=(nga,5,2,nt),dtype=float)
ntt=nt/12
rawdatasen3h=np.zeros(shape=(Nsen+1,nga,2,nz,ntt),dtype=float)
cloudpathsen3h=np.zeros(shape=(Nsen+1,nga,5,2,ntt),dtype=float)
#rawdata5003h=np.zeros(shape=(nga,2,nz,ntt),dtype=float)
#cloudpath5003h=np.zeros(shape=(nga,5,2,ntt),dtype=float)
alltitle=[]
allxdate=[]
allregion=[]
nt1=nz+nz+5+5
for iga in range(0,nga):
    name500m=CASE[iga]
    cdate=DATESTR[iga]
    cdatestr='%4i'%cdate[0]+'%2.2i'%cdate[1]+'%2.2i'%cdate[2]
    if name500m[0:3]=='MLY' :
        region=name500m[0:4]       
    else:
        region=name500m[0:3]          
    allregion.append(region)
    iskp=0
    fpath500=dirin500m+name500m+'_ALLTYPES_CLOUDWATER_PATH.TXT'
    ondedim500=readAscii(fpath, iskp)
    itt=0
    iks=0
    for it in range(0,nt):
        iks=it*nt1 
        print iks
        for iz in range(0,nz):
            ikk=iks+iz               
            rawdatasen[0,iga,0,iz,it]=ondedim500[ikk]               
            rawdatasen[0,iga,1,iz,it]=ondedim500[ikk+nz]
        iks=iks+2*nz         
        for itk in range(0,5):
            ikk=iks+itk
            cloudpathsen[0,iga,itk,0,it]=ondedim500[ikk]
            if ondedim500[ikk]<0:
                cloudpathsen[0,iga,itk,0,it]=0.
        iks=iks+5
        for itk in range(0,5):
            ikk=iks+itk
            cloudpathsen[0,iga,itk,1,it]=ondedim500[ikk]
            if ondedim500[ikk]<0:
                cloudpathsen[0,iga,itk,1,it]=0.
        for itk in range(0,2):
            for iz in range(0,nz):
                if rawdatasen[0,iga,itk,iz,it]<0:
                    rawdatasen[0,iga,0,iz,it]=None
        print iks
    for isn in range(1,Nsen+1):
        fpathsen=dirinsen+name3k+"%2.2i"%isn+'_ALLTYPES_CLOUDWATER_PATH.TXT'
        iskp=0
        ondedimsen=readAscii(fpath, iskp)
        itt=0
        iks=0
        for it in range(0,nt):
            iks=it*nt1 
            print iks
            for iz in range(0,nz):
                ikk=iks+iz               
                rawdatasen[isn,iga,0,iz,it]=ondedimsen[ikk]               
                rawdatasen[isn,iga,1,iz,it]=ondedimsen[ikk+nz]
            iks=iks+2*nz         
            for itk in range(0,5):
                ikk=iks+itk
                cloudpathsen[isn,iga,itk,0,it]=ondedimsen[ikk]
                if ondedimsen[ikk]<0:
                    cloudpathsen[isn,iga,itk,0,it]=0.
            iks=iks+5
            for itk in range(0,5):
                ikk=iks+itk
                cloudpathsen[isn,iga,itk,1,it]=ondedimsen[ikk]
                if ondedimsen[ikk]<0:
                    cloudpathsen[isn,iga,itk,1,it]=0.
            for itk in range(0,2):
                for iz in range(0,nz):
                    if rawdatasen[isn,iga,itk,iz,it]<0:
                        rawdatasen[isn,iga,0,iz,it]=None
            print iks    
    for itt in range(0,ntt):
        nnt1=itt*12-6
        nnt2=itt*12+6
        if itt==0:
            nnt1=0
            nnt2=7
        if itt==ntt-1:
            nnt1=nt-6
            nnt2=nt
        dnn=float(nnt2-nnt1)
        for iz in range(0,nz):
            for ik in range(0,2):                
                for isn in range(0,Nsen+1): 
                    rawdatasen3h[isn,iga,ik,iz,itt]=0.
                    for it1 in range(nnt1,nnt2):                                        
                        rawdatasen3h[isn,iga,ik,iz,itt]=rawdatasen3h[isn,iga,ik,iz,itt] \
                            +rawdatasen[isn,iga,ik,iz,it1]/dnn                   
        for itk in range(0,5):
            for ik in range(0,2):
                for isn in range(0,Nsen+1): 
                    cloudpathsen3h[isn,iga,itk,ik,itt]=0.
                    for it1 in range(nnt1,nnt2):                   
                        cloudpathsen3h[isn,iga,itk,ik,itt]=cloudpathsen3h[isn,iga,itk,ik,itt] \
                            +cloudpathsen[isn,iga,itk,ik,it1]/dnn
    alltitle.append(astr[iga]+' '+region)
    year=cdate[0]
    month=cdate[1]
    day=cdate[2]
    datestart=datetime.datetime(year,month,day,0,0,0)
    det=datetime.timedelta(hours=3)            
    dateiso=[]            
    for dt in range(0,ntt):
        dateiso.append(datestart+dt*det)
    xdate=[]              
    for tm in dateiso:
        xdate.append(datetime.datetime.strftime(tm,"%d/%b")) 
    allxdate.append(xdate)
###############################################################################
xdat=range(0,ntt)
fig,ax=plt.subplots(nrows=3,ncols=2,figsize=(18,12))
color_cycle=['deeppink', 'lime', 'b', 'y','indigo', 'cyan']
wd=[2,2,2,2,2]
jc=0
jr=0
ij=1
mycmap=cm.seismic
for iga in range(0,nga):
    region=allregion[iga]
    xdate=allxdate[iga]




    if jc==2:
        jc=0
        jr=jr+1
    mker=alltitle[iga]
    xdate=allxdate[iga]
    levs1=[0.2,0.4,0.6,0.8,1.0,1.2,1.5,2]
    levs2=[0.2,0.4,0.6,0.8,1.0,1.2,1.5,2]
    levs1=[0.0,0.01,0.03,0.04,0.06,0.09,0.12,0.15,0.18,0.21,0.24,0.27,0.3]
    levs2=[0.0,0.01,0.03,0.04,0.06,0.09,0.12,0.15,0.18,0.21,0.24,0.27,0.3]
    font = {'family' : 'serif',
        'color'  : 'k',
        'weight' : 'bold',
        'size'   : 16,
        }     
    plt.subplot(3,2,ij)
    zdat=rawdata5003h[iga,1,:,:]-rawdatasen3h[isn,iga,1,:,:]
    ax[jr,jc]=plt.contourf(xdat,ydat,zdat,cmap=mycmap ,extend='both', #colors='grey',
        linewidths=1.5,levels=levs2)  
#    ax[jr,jc]=plt.contourf(xdat,ydat,zdat, cmap=cm.jet,extend='both',
#        levels=levs2)#,linestyles=linetype2)  #linewidths=1.5,
    plt.axis([0, nt+1, 0, 16])
#    plt.clabel(ax[jr,jc],inline=1,fmt='%2.1f',fontsize=12)
    axx=fig.add_subplot(3,2,ij)
    ymajorLocator   = MultipleLocator(4)
    axx.yaxis.set_major_locator(ymajorLocator) 
    text1=mker  #r"($a$)"
    axx.text(1.5,16.5,text1,fontsize=18)                        
    axx.set_xticks(range(0,nt,96*4))
    xticklabels = [xdate[nn] for nn in range(0,nt,96*4)] 
    axx.set_xticklabels(xticklabels, size=16)
    plt.ylabel(r'Height ($km$)', fontdict=font)
    jc=jc+1
    ij=ij+1                
plt.show()
fig.subplots_adjust(left=0.1,bottom=0.1,right=1-0.1,top=1-0.1,hspace=0.4)
cax = fig.add_axes([0.2, 0.03, 0.6, 0.03])
fig.colorbar(ax[0,0], cax,extend='both',
              spacing='uniform', orientation='horizontal')
plt.savefig(dirpic+"rsolu_diff_ALLCASE_IceCloudWater.png",dpi=300)          
plt.show()
plt.close()     
                
                