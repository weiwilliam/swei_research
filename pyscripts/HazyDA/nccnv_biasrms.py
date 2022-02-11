#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 18 12:41:00 2018

@author: weiwilliam
"""
import os
import sys
sys.path.append('/Users/weiwilliam/AlbanyWork/Utility/Python3/functions')
import setuparea as setarea
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.dates import (DAILY, DateFormatter,
                              rrulewrapper, RRuleLocator)
import numpy as np
from datetime import datetime
from datetime import timedelta
import pandas as pd
import xarray as xa

mpl.rc('axes', titlesize=12,labelsize=12)
mpl.rc('xtick',labelsize=12)
mpl.rc('ytick',labelsize=12)
mpl.rc('legend',fontsize='large')
minussign=u'\u2212'
#mpl.rc('lines',linewidth=1.2)

outpath='/Users/weiwilliam/AlbanyWork/R2O_Project/images/tmp/Conventional/BiasRMS'
path='/Volumes/WD2TB/ResearchData/2_R2O/nc_DiagFiles'

varlist=['t'] #['ps','sst','gps','q','t','uv','tcp']
unitlist=['K'] #['mb','K','%','g/kg','K','m/s','mb']
explist=np.array(['prctrl','praero'])
expnlist=['CTL_cyc','AER_cyc']
enum=explist.shape[0]

sdate=2017080100
edate=2017082818
hint=6

loop='anl' #ges,anl
if (loop=='ges'):
   lpstr='OMF'
elif (loop=='anl'):
   lpstr='OMA'

area='TRO'
minlon, maxlon, minlat, maxlat, crosszero, cyclic=setarea.setarea(area)
print(area,minlat,maxlat,minlon,maxlon,crosszero)
fsave=1

useqc=1
if (useqc):
    qcflg='qc'
else:
    qcflg='noqc'

zpltlst=[0]

imgsavpath=outpath+'/'+expnlist[1]+'/'+area
if ( not os.path.exists(imgsavpath) ):
   os.makedirs(imgsavpath)

def ndate(cdate,hinc):
    yy=int(str(cdate)[:4])
    mm=int(str(cdate)[4:6])
    dd=int(str(cdate)[6:8])
    hh=int(str(cdate)[8:10])
    dstart=datetime(yy,mm,dd,hh)
    dnew=dstart+timedelta(hours=hinc)
    dnewint=int(str('%4.4d' % dnew.year)+str('%2.2d' % dnew.month)+
                str('%2.2d' %dnew.day)+str('%2.2d' % dnew.hour))
    return dnewint

syy=int(str(sdate)[:4]); smm=int(str(sdate)[4:6])
sdd=int(str(sdate)[6:8]); shh=int(str(sdate)[8:10])
eyy=int(str(edate)[:4]); emm=int(str(edate)[4:6])
edd=int(str(edate)[6:8]); ehh=int(str(edate)[8:10])

date1 = datetime(syy,smm,sdd,shh)
date2 = datetime(eyy,emm,edd,ehh)
delta = timedelta(hours=6)
dates = pd.date_range(start=date1, end=date2, freq=delta)

xdate2= date2+delta
xdates= mdates.drange(date1, xdate2, delta)

rule = rrulewrapper(DAILY, byhour=6, interval=5)
loc = RRuleLocator(rule)
formatter = DateFormatter('%Y%h %n %d %Hz')

# Calculate how many cases are going to do statistic.
tnum=0
dlist=[]
cdate=sdate
while (cdate<=edate):
    dlist.append(str(cdate))
    tnum=tnum+1
    cdate=ndate(cdate,6)
    
print('Total cases number is %d' % tnum )

ptop=np.array((1000.,900.,800.,600.,400.,300.,250.,200.,150.,100.,50.,0.))
pbot=np.array((1200.,1000.,900.,800.,600.,400.,300.,250.,200.,150.,100.,50.))
znum=ptop.size

uidx=0
for var in varlist:
    unit=unitlist[uidx]
    if (var=='ps' or var=='sst' or var=='tcp'):
        icount=np.zeros((tnum,2))
    else:
        icount=np.zeros((tnum,znum,2))
    d=0
    for date in dlist:
        cnvdfile='diag_conv_'+var+'_'+loop+'.'+date+'.nc'
        print('File: %s' % cnvdfile)
        infile1=path+'/'+explist[0]+'/'+date+'/'+cnvdfile
        infile2=path+'/'+explist[1]+'/'+date+'/'+cnvdfile
        if (os.path.exists(infile1) and os.path.exists(infile2)):
            print('Processing Cnvfile: %s' %(cnvdfile))
            ds1=xa.open_dataset(infile1)
            ds2=xa.open_dataset(infile2)
            #print('Load Data Elapsed: %f [s]' %(end-start))
            try:
                omg_mean
            except NameError:
                if (var=='ps' or var=='sst' or var=='tcp'):
                    omg_mean=np.zeros((tnum,2),dtype='float')
                    omg_rmsq=np.zeros_like(omg_mean)
                    omg_mean[:,:]=np.nan
                    omg_rmsq[:,:]=np.nan
                else:
                    omg_mean=np.zeros((tnum,znum,2),dtype='float')
                    omg_rmsq=np.zeros_like(omg_mean)
                    omg_mean[:,:,:]=np.nan
                    omg_rmsq[:,:,:]=np.nan
        else:
            print('%s is not existing'%(cnvdfile))
            d=d+1
            continue
        
        rlat1=ds1.obsdata[2,:]
        rlon1=ds1.obsdata[3,:]
        rlat2=ds2.obsdata[2,:]
        rlon2=ds2.obsdata[3,:]
        
        if (crosszero):
            rlon1[rlon1>=maxlon]=rlon1[rlon1>=maxlon]-360.
            rlon2[rlon2>=maxlon]=rlon2[rlon2>=maxlon]-360.
        
        pres1=ds1.obsdata[5,:]
        iuse1=ds1.obsdata[11,:]
        pres2=ds2.obsdata[5,:]
        iuse2=ds2.obsdata[11,:]
        
        mask1=(ds1.vartype==var)
        mask2=(ds2.vartype==var)
        
        if (area!='Glb'):
            mask1=(mask1)&((rlon1<maxlon)&(rlon1>minlon)&(rlat1>minlat)&(rlat1<maxlat))
            mask2=(mask2)&((rlon2<maxlon)&(rlon2>minlon)&(rlat2>minlat)&(rlat2<maxlat))
        
        if (useqc):
            mask1=(mask1)&(iuse1==1)
            mask2=(mask2)&(iuse2==1)
        
        if (var=='ps' or var=='sst' or var=='tcp'):
            dpar1=ds1.obsdata[17,:]
            dpar2=ds2.obsdata[17,:]
            
            icount[d,0]=np.count_nonzero(mask1)
            icount[d,1]=np.count_nonzero(mask2)
            
            dpar1=xa.where(mask1,dpar1,np.nan)
            dpar2=xa.where(mask2,dpar2,np.nan)
            
            omg_mean[d,0]=np.nanmean(dpar1)
            omg_mean[d,1]=np.nanmean(dpar2)
            
            omg_rmsq[d,0]=np.sqrt(np.nanmean(np.square(dpar1)))
            omg_rmsq[d,1]=np.sqrt(np.nanmean(np.square(dpar2)))
            
        else:
            if (var=='uv'):
                dpar1=ds1.obsdata[17,:]
                dpar2=ds2.obsdata[17,:]
            else:
                dpar1=ds1.obsdata[17,:]
                dpar2=ds2.obsdata[17,:]
                
            for z in np.arange(znum):
                zmask1=(mask1)&((pres1<pbot[z])&(pres1>ptop[z]))
                zmask2=(mask2)&((pres2<pbot[z])&(pres2>ptop[z]))
                
                icount[d,z,0]=np.count_nonzero(zmask1)
                icount[d,z,1]=np.count_nonzero(zmask2)
                
                zdpar1=xa.where(zmask1,dpar1,np.nan)
                zdpar2=xa.where(zmask2,dpar2,np.nan)
                
                omg_mean[d,z,0]=np.nanmean(zdpar1)
                omg_mean[d,z,1]=np.nanmean(zdpar2)
            
                omg_rmsq[d,z,0]=np.sqrt(np.nanmean(np.square(zdpar1)))
                omg_rmsq[d,z,1]=np.sqrt(np.nanmean(np.square(zdpar2)))
        d=d+1
                
    if (var=='ps' or var=='sst' or var=='tcp'):
        fig,ax=plt.subplots(2,1,sharex=True,figsize=(9,3.8))
        fig.subplots_adjust(hspace=0)
        for a in np.arange(2):
            ax[a].set_prop_cycle(color=['red','blue'])
            ax[a].grid()
        
        ax[1].plot_date(xdates,omg_mean,'-')
        ax[0].xaxis.set_major_locator(loc)
        ax[0].xaxis.set_major_formatter(formatter)
        ax[0].xaxis.set_tick_params(rotation=30, labelsize=10)
        ax[0].set_title('%s %s[%s]' %(area,var.upper(),unit))
        ax[0].plot_date(xdates,omg_rmsq,'--')
        ax[0].set_ylabel('RMS %s [%s]' %(lpstr,unitlist[uidx]))
        ax[1].set_ylabel('Mean %s [%s]'%(lpstr,unitlist[uidx]))
        lglist=np.zeros((2,2),dtype='<U30')
        for ex in np.arange(2):
            lglist[0,ex]=expnlist[ex]+'(%.2f)' %(np.nanmean(omg_rmsq[:,ex]))
            lglist[1,ex]=expnlist[ex]+'(%.2f)' %(np.nanmean(omg_mean[:,ex]))
        ax[0].legend(lglist[0,:])
        ax[1].legend(lglist[1,:])
        if (fsave):
            fig.savefig(imgsavpath+'/%s_%s_%s_%s_%s_%s_BIASRMS.png' 
                        %(area,loop,var,explist[0],explist[1],qcflg), dpi=200)
            plt.close()
        
        fig,ax=plt.subplots(figsize=(9,3.8))
        fig.subplots_adjust(hspace=0)
        ax.set_prop_cycle(color=['red','blue'])
        ax.grid()
        ax.xaxis.set_major_locator(loc)
        ax.xaxis.set_major_formatter(formatter)
        ax.xaxis.set_tick_params(rotation=30, labelsize=10)
        ax.set_title('%s %s' %(area,var.upper()))
        ax.plot_date(xdates,icount,'-')
        ax.set_ylabel('OBS number')
        lglist=np.zeros((2),dtype='<U30')
        for ex in np.arange(2):
            lglist[ex]=expnlist[ex]+'(%8.4f)' %(np.nanmean(icount[:,ex]))
        ax.legend(lglist[:])
        if (fsave):
            fig.savefig(imgsavpath+'/%s_%s_%s_%s_%s_%s_OBSnum.png' 
                        %(area,loop,var,explist[0],explist[1],qcflg), dpi=200)
            plt.close()
        
    else:
        fig,ax=plt.subplots(1,3,sharey=True)
        fig.subplots_adjust(wspace=0.05)
        ax[0].set_prop_cycle(color=['red','blue'])
        ax[1].set_prop_cycle(color=['red','blue'])
        ax[2].set_prop_cycle(color=['k','k'],linestyle=['-','--'])
        biasplot=np.nanmean(omg_mean,axis=0)
        rmsplot=np.nanmean(omg_rmsq,axis=0)
        diffplot=np.zeros_like(biasplot)
        diffplot[:,0]=np.diff(biasplot,axis=1)[:,0]
        diffplot[:,1]=np.diff(rmsplot,axis=1)[:,0]
        ax[0].plot(biasplot,pbot,'-')
        ax[0].set_title('Mean %s [%s]'%(lpstr,unitlist[uidx]))
        ax[1].plot(rmsplot,pbot,'--')
        ax[1].set_title('RMS %s [%s]'%(lpstr,unitlist[uidx]))
        ax[2].plot(diffplot,pbot)
        ax[2].set_title(expnlist[1]+'-'+expnlist[0])
        ax[0].invert_yaxis()
        fig.suptitle('%s [%s]' %(var.upper(),unit))
        ax[0].legend(expnlist)
        ax[2].legend(['Mean %s'%(lpstr),'RMS %s'%(lpstr)])
        ax[0].grid()
        ax[1].grid()
        ax[2].grid()
        if (fsave):
            fig.savefig(imgsavpath+'/%s_%s_%s_%s_%s_%s_BIASRMS.png' 
                        %(area,loop,var,explist[0],explist[1],qcflg), dpi=200)
            plt.close()
        
        for z in zpltlst:
            fig,ax=plt.subplots(2,1,sharex=True,figsize=(9,3.8))
            fig.subplots_adjust(hspace=0)
            for a in np.arange(2):
                ax[a].set_prop_cycle(color=['red','blue'])
                ax[a].grid()
        
            ax[1].plot_date(dates,omg_mean[:,z,:],'-')
            ax[0].xaxis.set_major_locator(loc)
            ax[0].xaxis.set_major_formatter(formatter)
            ax[0].xaxis.set_tick_params(rotation=30, labelsize=10)
            ax[0].set_title('%s %s [%s] %.1f %s %.1f [hPa]' %(area,var.upper(),unit,ptop[z],minussign,pbot[z]),loc='left')
            ax[0].plot_date(dates,omg_rmsq[:,z,:],'--')
            ax[0].set_ylabel('RMS %s [%s]'%(lpstr,unitlist[uidx]))
            ax[1].set_ylabel('Mean %s [%s]'%(lpstr,unitlist[uidx]))
            lglist=np.zeros((2,2),dtype='<U30')
            for ex in np.arange(2):
                lglist[0,ex]=expnlist[ex]+'(%.2f)' %(np.nanmean(omg_rmsq[:,z,ex]))
                lglist[1,ex]=expnlist[ex]+'(%.2f)' %(np.nanmean(omg_mean[:,z,ex]))
            ax[0].legend(lglist[0,:])
            ax[1].legend(lglist[1,:])
            
            if (fsave):
                fig.savefig(imgsavpath+'/%s_%s_%s_%s_%s_%i_%s_BIASRMS.png' 
                            %(area,loop,var,explist[0],explist[1],pbot[z],qcflg), dpi=200)
                plt.close()
                
            fig,ax=plt.subplots(figsize=(9,3.8),constrained_layout=1)
            #fig.subplots_adjust(hspace=0)
            ax.set_prop_cycle(color=['red','blue'])
            ax.grid()
            ax.xaxis.set_major_locator(loc)
            ax.xaxis.set_major_formatter(formatter)
            ax.xaxis.set_tick_params(rotation=30, labelsize=10)
            ax.set_title('%s %s %.1f - %.1f' %(area,var.upper(),pbot[z],ptop[z]))
            ax.plot_date(xdates,icount[:,z,:],'-')
            ax.set_ylabel('OBS number')
            lglist=np.zeros((2),dtype='<U30')
            for ex in np.arange(2):
                lglist[ex]=expnlist[ex]+'(%8.2f)' %(np.nanmean(icount[:,z,ex]))
            ax.legend(lglist[:])
            if (fsave):
                fig.savefig(imgsavpath+'/%s_%s_%s_%s_%s_P%i_%s_OBSnum.png'
                            %(area,loop,var,explist[0],explist[1],pbot[z],qcflg), dpi=200)
                plt.close()
                
    uidx=uidx+1
    

