#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 18 12:41:00 2018

@author: weiwilliam
"""
import sys, os, platform
import numpy as np
import xarray as xa
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mpcrs
import matplotlib.dates as mdates
from matplotlib.dates import (DAILY, DateFormatter,
                              rrulewrapper, RRuleLocator)
import setuparea as setarea
from plot_utils import setupax_2dmap, set_size
from utils import ndate,setup_cmap
from stats import calculate_stats, is_significant
from datetime import datetime, timedelta

tlsize=12 ; lbsize=10
mpl.rc('axes', titlesize=tlsize,labelsize=lbsize)
mpl.rc('xtick',labelsize=lbsize)
mpl.rc('ytick',labelsize=lbsize)
mpl.rc('legend',fontsize='small')
plot_ts = 0
fsave=1 ; ffmt='png' ; ptsize=4
axe_w=3 ; axe_h=3 ; quality=600
diagsuffix='nc4'
minussign=u'\u2212'
#mpl.rc('lines',linewidth=1.2)

outputpath = '/data/users/swei/ResearchData/HazySky_DA/Images/DiagFiles/conv/'
ncdiagpath = '/data/users/swei/ResearchData/HazySky_DA/AeroObsStats/nc_diag'

var_unit_mapping = {
    'ps': 'mb',
    'sst': 'K',
    'gps': '%',
    'q': 'g/kg',
    't': 'K',
    'uv': 'm/s',
    'tcp': 'mb',
}


plotvars=['t', 'q']
bufrtype='all' # SST: 181-199
explist=np.array(['hazyda_ctrl','hazyda_aero'])
expnlist=['CTL','AER']
nexp=explist.size
exps_fname_str=''
for exp in expnlist:
    exps_fname_str += exp+'_'
exps_fname_str=exps_fname_str[:-1]

sdate=2020061000
edate=2020071018
hint=6

loop='ges' #ges,anl
area='TRO'
useqc=1

if (loop=='ges'):
   lpstr='OMF'
elif (loop=='anl'):
   lpstr='OMA'
if (useqc):
    qcflg='qc'
else:
    qcflg='noqc'
minlon, maxlon, minlat, maxlat, crosszero, cyclic=setarea.setarea(area)
print(area,minlat,maxlat,minlon,maxlon,crosszero)

zpltlst=[0,1,2,3,4,5,6,7,8]

imgsavpath=outputpath+'/'+area
if ( not os.path.exists(imgsavpath) ):
   os.makedirs(imgsavpath)

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
formatter = DateFormatter('%Y %h %n %d %Hz')

# Calculate how many cases are going to do statistic.
tnum=0
dlist=[]
cdate=sdate
while (cdate<=edate):
    dlist.append(str(cdate))
    tnum=tnum+1
    cdate=ndate(hint,cdate)
    
print('Total cases number is %d' % tnum )

ptop=np.array((1000.,900.,800.,600.,400.,300.,250.,200.,150.,100.,50.,0.))
pbot=np.array((1200.,1000.,900.,800.,600.,400.,300.,250.,200.,150.,100.,50.))
pmid=0.5*(ptop+pbot)
znum=ptop.size
zlabels=[]
for z in np.arange(znum):
    zlabels.append('%.0i-%.0i'%(ptop[z],pbot[z]))

for var in plotvars:
    unit = var_unit_mapping[var]
    icount=np.zeros((tnum,znum,2))
    d=0
    for date in dlist:
        cnvdfile = 'diag_conv_'+var+'_'+loop+'.'+date+'.'+diagsuffix
        infile1 = ncdiagpath+'/'+explist[0]+'/'+date+'/'+cnvdfile
        infile2 = ncdiagpath+'/'+explist[1]+'/'+date+'/'+cnvdfile
        if (os.path.exists(infile1) and os.path.exists(infile2)):
            print('Processing Cnvfile: %s' %(cnvdfile))
            ds1=xa.open_dataset(infile1)
            ds2=xa.open_dataset(infile2)
            try:
                omg_mean
            except NameError:
                omg_mean=np.zeros((tnum, znum, 2), dtype=np.float32)
                omg_rmsq=np.zeros_like(omg_mean)
                t_signif = np.zeros((tnum, znum), dtype=np.int32)
                z_signif = np.zeros((znum), dtype=np.int32)
                omg_mean[:,:,:] = np.nan
                omg_rmsq[:,:,:] = np.nan
        else:
            print('%s is not existing'%(cnvdfile))
            d=d+1
            continue
        
        rlat1=ds1.Latitude.data
        rlat2=ds2.Latitude.data
        rlon1=(ds1.Longitude.data+180)%360-180
        rlon2=(ds2.Longitude.data+180)%360-180
        
        pres1=ds1.Pressure
        pres2=ds2.Pressure

        iuse1=ds1.Analysis_Use_Flag
        iuse2=ds2.Analysis_Use_Flag
        type1=ds1.Observation_Type
        type2=ds2.Observation_Type
        
        mask1=(~np.isnan(ds1.nobs))
        mask2=(~np.isnan(ds2.nobs))
        
        if (area!='Glb'):
            mask1=(mask1)&((rlon1<maxlon)&(rlon1>minlon)&(rlat1>minlat)&(rlat1<maxlat))
            mask2=(mask2)&((rlon2<maxlon)&(rlon2>minlon)&(rlat2>minlat)&(rlat2<maxlat))
        
        if (useqc):
            mask1=(mask1)&(iuse1==1)
            mask2=(mask2)&(iuse2==1)

        if (bufrtype!='all'):
            mask1=(mask1)&(type1==int(bufrtype))
            mask2=(mask2)&(type2==int(bufrtype))
        
        if (var=='uv'):
            dpar1=ds1.u_Obs_Minus_Forecast_adjusted
            dpar2=ds2.u_Obs_Minus_Forecast_adjusted
        else:
            dpar1=ds1.Obs_Minus_Forecast_adjusted
            dpar2=ds2.Obs_Minus_Forecast_adjusted
            
        for z in np.arange(znum):
            zmask1=(mask1)&((pres1<pbot[z])&(pres1>ptop[z]))
            zmask2=(mask2)&((pres2<pbot[z])&(pres2>ptop[z]))
            
            icount[d,z,0]=np.count_nonzero(zmask1)
            icount[d,z,1]=np.count_nonzero(zmask2)
            
            zdpar1 = dpar1[zmask1]
            zdpar2 = dpar2[zmask2]
            
            omg_mean[d,z,0]=np.nanmean(zdpar1)
            omg_mean[d,z,1]=np.nanmean(zdpar2)
        
            omg_rmsq[d,z,0]=np.sqrt(np.nanmean(np.square(zdpar1)))
            omg_rmsq[d,z,1]=np.sqrt(np.nanmean(np.square(zdpar2)))
            t_signif[d, z] = is_significant(zdpar1, zdpar2, 10000)
        d+=1
    
    fig,ax=plt.subplots(1,3,sharey=True)
    fig.subplots_adjust(left=0.2,right=0.95,wspace=0.05)
    ax[0].set_prop_cycle(color=['blue','red'],linestyle=['-','-'],marker=['o','o'])
    ax[1].set_prop_cycle(color=['blue','red'],linestyle=['--','--'],marker=['o','o'])
    ax[2].set_prop_cycle(color=['k','k'],linestyle=['-','--'],marker=['o','o'])

    biasplot=np.nanmean(omg_mean,axis=0)
    rmsplot=np.nanmean(omg_rmsq,axis=0)
    diffplot=np.zeros_like(biasplot)
    diffplot[:,0]=np.diff(biasplot,axis=1)[:,0]
    diffplot[:,1]=np.diff(rmsplot,axis=1)[:,0]

    for z in np.arange(znum):
        z_signif[z] = is_significant(omg_rmsq[:, z, 0], omg_rmsq[:, z, 1], 10000)
    z_sigmask = z_signif!=0
    z_sigplotval = np.full(z_signif[z_sigmask].shape, diffplot[:, 1].min()-abs(diffplot[:, 1].min())*0.1)

    ax[0].plot(biasplot,pmid,lw=1,ms=1.5)
    ax[0].set_title('Mean %s [%s]'%(lpstr,unit))
    ax[1].plot(rmsplot,pmid,lw=1,ms=1.5)
    ax[1].set_title('RMS %s [%s]'%(lpstr,unit))
    ax[2].plot(diffplot,pmid,lw=1,ms=1.5)
    ax[2].plot(z_sigplotval, pmid[z_sigmask], 's', c='g', ms=1.)
    ax[2].set_title('%s%s%s [%s]'%(expnlist[1],minussign,expnlist[0],unit))
    ax[0].invert_yaxis()
    ax[0].set_yscale('log')
    ax[0].set_yticks(pmid[::2])
    ax[0].set_yticklabels(zlabels[::2])
    ax[0].set_ylabel('Pressure [hPa]')
    fig.suptitle('%s [%s]' %(var.upper(),unit))
    ax[0].legend(expnlist)
    ax[2].legend(['Mean','RMS'],loc=2)
    ax[0].grid()
    ax[1].grid()
    ax[2].grid()
    if (fsave):
        fname=('%s/%s_%s_%s_%s_%s_bufr%s_BIASRMS.%s_%s.png'
               %(imgsavpath,area,loop,var,exps_fname_str,qcflg,bufrtype,sdate,edate))
        print(fname,flush=1)
        fig.savefig(fname, dpi=quality)
        plt.close()
    
    if plot_ts:
        for z in zpltlst:
            fig,ax=plt.subplots(2,1,sharex=True,figsize=(9,3.8))
            fig.subplots_adjust(hspace=0.1)
            for a in np.arange(2):
               ax[a].set_prop_cycle(color=['blue','red'])
               ax[a].grid()
        
            ax[1].plot_date(dates,omg_mean[:,z,:],'-o',ms=1.5,lw=1)
            ax[0].xaxis.set_major_locator(loc)
            ax[0].xaxis.set_major_formatter(formatter)
            ax[0].xaxis.set_tick_params(labelsize=10)
            ax[0].set_title('%s %s [%s] %.1f %s %.1f [hPa]' %(area,var.upper(),unit,ptop[z],minussign,pbot[z]),loc='left')
            ax[0].plot_date(dates,omg_rmsq[:,z,:],'--o',ms=1.5,lw=1)
            ax[0].set_ylabel('RMS %s [%s]'%(lpstr,unit))
            ax[1].set_ylabel('Mean %s [%s]'%(lpstr,unit))
            lglist=np.zeros((2,2),dtype='<U30')
            for ex in np.arange(nexp):
               lglist[0,ex]=expnlist[ex]+'(%.2f)' %(np.nanmean(omg_rmsq[:,z,ex]))
               lglist[1,ex]=expnlist[ex]+'(%.2f)' %(np.nanmean(omg_mean[:,z,ex]))
            ax[0].legend(lglist[0,:])
            ax[1].legend(lglist[1,:])
            
            if (fsave):
               fname=('%s/%s_%s_%s_%s_%i_%s_bufr%s_BIASRMS.%s_%s.png' 
                      %(imgsavpath,area,loop,var,exps_fname_str,pbot[z],qcflg,bufrtype,sdate,edate))
               print(fname,flush=1)
               fig.savefig(fname, dpi=quality)
               plt.close()
                
