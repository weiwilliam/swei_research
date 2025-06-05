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
plot = 1
fsave=1 ; ffmt='png' ; ptsize=4
axe_w=3 ; axe_h=3 ; quality=600
diagsuffix='nc4'
minussign=u'\u2212'
#mpl.rc('lines',linewidth=1.2)

outputpath = '/data/users/swei/ResearchData/HazySky_DA/Images/DiagFiles/conv/'
ncdiagpath = '/data/users/swei/ResearchData/HazySky_DA/AeroObsStats/nc_diag'

sfctype_list=['180','181','182','183','187']
varlist=['sst'] #['ps','sst','gps','q','t','uv','tcp']
unitlist=['K'] #['mb','K','%','g/kg','K','m/s','mb']
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
useqc=0

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

uidx=0
for var in varlist:
    unit=unitlist[uidx]
    icount=np.zeros((tnum,2))
    d=0
    for date in dlist:
        cnvdfile='diag_conv_'+var+'_'+loop+'.'+date+'.'+diagsuffix
        infile1=ncdiagpath+'/'+explist[0]+'/'+date+'/'+cnvdfile
        infile2=ncdiagpath+'/'+explist[1]+'/'+date+'/'+cnvdfile
        if (os.path.exists(infile1) and os.path.exists(infile2)):
            print('Processing Cnvfile: %s' %(cnvdfile))
            ds1=xa.open_dataset(infile1)
            ds2=xa.open_dataset(infile2)
            try:
                omg_mean
            except NameError:
                omg_mean = np.zeros((tnum, 2), dtype=np.float32)
                omg_rmsq = np.zeros_like(omg_mean)
                corr = np.zeros_like(omg_mean)
                significant = np.zeros((tnum), dtype=np.int32)
                omg_mean[:, :] = np.nan
                omg_rmsq[:, :] = np.nan
                corr[:] = np.nan
        else:
            print('%s is not existing'%(cnvdfile))
            d += 1
            continue
        
        rlat1 = ds1.Latitude.data
        rlat2 = ds2.Latitude.data
        rlon1 = (ds1.Longitude.data+180)%360-180
        rlon2 = (ds2.Longitude.data+180)%360-180
        
        iuse1 = ds1.Analysis_Use_Flag
        iuse2 = ds2.Analysis_Use_Flag
        type1 = ds1.Observation_Type
        type2 = ds2.Observation_Type
        
        mask1 = (~np.isnan(ds1.nobs))
        mask2 = (~np.isnan(ds2.nobs))
        
        if (area!='Glb'):
            mask1=(mask1)&((rlon1<maxlon)&(rlon1>minlon)&(rlat1>minlat)&(rlat1<maxlat))
            mask2=(mask2)&((rlon2<maxlon)&(rlon2>minlon)&(rlat2>minlat)&(rlat2<maxlat))
        
        if (useqc):
            mask1=(mask1)&(iuse1==1)
            mask2=(mask2)&(iuse2==1)

        if (bufrtype!='all'):
            mask1=(mask1)&(type1==int(bufrtype))
            mask2=(mask2)&(type2==int(bufrtype))
        
        obs1 = ds1.Observation.data[mask1]
        obs2 = ds2.Observation.data[mask2]
        hfx1 = obs1 - ds1.Obs_Minus_Forecast_adjusted.data[mask1]
        hfx2 = obs2 - ds2.Obs_Minus_Forecast_adjusted.data[mask2]
        dpar1 = ds1.Obs_Minus_Forecast_adjusted.data[mask1]
        dpar2 = ds2.Obs_Minus_Forecast_adjusted.data[mask2]
        
        icount[d, 0] = np.count_nonzero(mask1)
        icount[d, 1] = np.count_nonzero(mask2)
        
        omg_mean[d, 0] = np.nanmean(dpar1)
        omg_mean[d, 1] = np.nanmean(dpar2)
        
        omg_rmsq[d, 0] = np.sqrt(np.nanmean(np.square(dpar1)))
        omg_rmsq[d, 1] = np.sqrt(np.nanmean(np.square(dpar2)))

        significant[d] = is_significant(dpar1, dpar2, 10000)
        stats_dict1 = calculate_stats(obs1, hfx1)
        stats_dict2 = calculate_stats(obs2, hfx2)
        corr[d, 0] = stats_dict1['R']
        corr[d, 1] = stats_dict2['R']

        if d==0:
            all_obs1 = obs1
            all_obs2 = obs2
            all_hfx1 = hfx1
            all_hfx2 = hfx2
        else:
            all_obs1 = np.concatenate((all_obs1, obs1))
            all_obs2 = np.concatenate((all_obs2, obs2))
            all_hfx1 = np.concatenate((all_hfx1, hfx1))
            all_hfx2 = np.concatenate((all_hfx2, hfx2))
            
        d+=1

    print('Significant:', significant)
    sigmask = significant!=0
    sigplotval = np.full(xdates[sigmask].shape, omg_rmsq.min()*0.9)
    print(explist[0])
    print('Correlation:', corr[:, 0])
    calculate_stats(all_obs1, all_hfx1, verbose=True)
    print(explist[1])
    print('Correlation:', corr[:, 1])
    calculate_stats(all_obs2, all_hfx2, verbose=True)

    if plot:
        fig, ax = plt.subplots(2,1,sharex=True,figsize=(9,3.8))
        fig.subplots_adjust(hspace=0.1)
        for a in np.arange(2):
            ax[a].set_prop_cycle(color=['blue','red','green'])
            ax[a].grid()
        
        ax[1].plot_date(xdates,omg_mean,'-o',lw=1,ms=1.5)
        ax[0].xaxis.set_major_locator(loc)
        ax[0].xaxis.set_major_formatter(formatter)
        ax[0].xaxis.set_tick_params(labelsize=10)
        ax[0].set_title('%s %s[%s]' %(area,var.upper(),unit),loc='left')
        ax[0].plot_date(xdates, omg_rmsq, '--o', lw=1, ms=1.5)
        ax[0].plot_date(xdates[sigmask], sigplotval, 's', ms=1.)
        ax[0].set_ylabel('RMS %s [%s]' %(lpstr,unitlist[uidx]))
        ax[1].set_ylabel('Mean %s [%s]'%(lpstr,unitlist[uidx]))
        lglist=np.zeros((2,2),dtype='<U30')
        for ex in np.arange(nexp):
            lglist[0,ex] = expnlist[ex]+'(%.2f)' %(np.nanmean(omg_rmsq[:,ex]))
            lglist[1,ex] = expnlist[ex]+'(%.2f)' %(np.nanmean(omg_mean[:,ex]))
        ax[0].legend(lglist[0,:])
        ax[1].legend(lglist[1,:])
        if (fsave):
            fname='%s/%s_%s_%s_%s_%s_bufr%s_BIASRMS.%s_%s.png'%(imgsavpath,area,loop,var,exps_fname_str,qcflg,bufrtype,sdate,edate)
            print(fname,flush=1)
            fig.savefig(fname, dpi=quality)
            plt.close()
        
    uidx=uidx+1
    

