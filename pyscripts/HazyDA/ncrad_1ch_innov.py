#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 21:59:35 2019

@author: weiwilliam
"""
import os, sys, platform
os_name=platform.system()
if (os_name=='Darwin'):
    rootpath='/Users/weiwilliam'
    rootarch='/Volumes/WD2TB/ResearchData'
elif (os_name=='Windows'):
    rootpath='F:\GoogleDrive_NCU\Albany'
    rootarch='F:\ResearchData'
    rootgit='F:\GitHub\swei_research'
elif (os_name=='Linux'):
    if (os.path.exists('/scratch1')):
        rootpath='/scratch2/BMC/gsd-fv3-dev/Shih-wei.Wei'
        rootarch='/scratch2/BMC/gsd-fv3-dev/Shih-wei.Wei/ResearchData'
        rootgit='/home/Shih-wei.Wei/research'
    elif (os.path.exists('/glade')):
        rootpath='/glade/work/swei/output/images'
        rootarch='/scratch2/BMC/gsd-fv3-dev/Shih-wei.Wei/ResearchData'
        rootgit='/glade/u/home/swei/research'
        machine='Cheyenne'
    elif (os.path.exists('/cardinal')):
        rootpath='/data/users/swei/Images'
        rootarch='/scratch/users/swei/ncdiag'
        rootgit='/home/swei/research'
        machine='S4'
sys.path.append(rootgit+'/pyscripts/functions')
import setuparea as setarea
from utils import ndate 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mpcrs
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from datetime import datetime
from datetime import timedelta
import matplotlib as mpl
import matplotlib.dates as mdates
from matplotlib.dates import (DAILY, DateFormatter,
                              rrulewrapper, RRuleLocator)
import matplotlib.ticker as ticker
import xarray as xa
import seaborn as sb
import pandas as pd

tlsize=8 ; lbsize=8
mpl.rc('axes', titlesize=tlsize,labelsize=lbsize)
mpl.rc('xtick',labelsize=lbsize)
mpl.rc('ytick',labelsize=lbsize)
mpl.rc('legend',fontsize='small')
fsave=1 ; ffmt='png' ; ptsize=4
axe_w=6 ; axe_h=2.7 ; quality=300
diagsuffix='nc4'

outputpath=rootpath+'/HazyDA'
inputpath=rootarch

explist=np.array(['hazyda_ctrl','hazyda_aero'])
expnlist=['CTL','AER']
sensor='iasi_metop-a'

chkwvl=13.
colormap=['bwr','bwr','bwr']

sdate=2020061000
edate=2020061000
hint=6
spectral_range=slice(600,1300)

area='r2o10'
minlon, maxlon, minlat, maxlat, crosszero, cyclic=setarea.setarea(area)
print(area,minlat,maxlat,minlon,maxlon,crosszero,cyclic)

loop='ges' #ges,anl
if loop=='anl':
    tlstr='OMA'
elif loop=='ges':
    tlstr='OMF'
usebc=0
if (usebc):
    bcflg='bc'
else:
    bcflg='nobc'
useqc=1
if (useqc):
    qcflg='qc'
else:
    qcflg='noqc'
    
rmflier=1
wateronly=0
if (wateronly):
    waterflg='water'
else:
    waterflg='all'

imgsavpath=outputpath+'/'+expnlist[1]+'/'+area
if ( not os.path.exists(imgsavpath) ):
   os.makedirs(imgsavpath)
    
syy=int(str(sdate)[:4]); smm=int(str(sdate)[4:6])
sdd=int(str(sdate)[6:8]); shh=int(str(sdate)[8:10])
eyy=int(str(edate)[:4]); emm=int(str(edate)[4:6])
edd=int(str(edate)[6:8]); ehh=int(str(edate)[8:10])

date1 = datetime(syy,smm,sdd,shh)
date2 = datetime(eyy,emm,edd,ehh)
delta = timedelta(hours=hint)

dates = pd.date_range(start=date1, end=date2, freq=delta)

xdate2= date2+delta
xdates= mdates.drange(date1, xdate2, delta)

rule = rrulewrapper(DAILY, byhour=6, interval=5)
loc = RRuleLocator(rule)
formatter = DateFormatter('%Y%h %n %d %Hz')

tnum=0
dlist=[]
cdate=sdate
while (cdate<=edate):
    dlist.append(str(cdate))
    tnum=tnum+1
    cdate=ndate(hint,cdate)

dates_count=0
for date in dlist:
    raddfile='diag_'+sensor+'_'+loop+'.'+date+'.'+diagsuffix
    infile0=inputpath+'/'+explist[0]+'/'+date+'/'+raddfile
    infile1=inputpath+'/'+explist[1]+'/'+date+'/'+raddfile
    
    if (os.path.exists(infile0) and
        os.path.exists(infile1) ):
        print('Processing Radfile: %s' %(raddfile))
        ds0=xa.open_dataset(infile0)
        ds1=xa.open_dataset(infile1)
        # npts0=ds0.obsloc.size
        npts0=int(ds0.nobs.size/ds0.nchans.size)
        npts1=int(ds1.nobs.size/ds1.nchans.size)
        nchs0=ds0.nchans.size
        nchs1=ds1.nchans.size
        ds0=ds0.swap_dims({"nchans":"wavenumber"})
        ds1=ds1.swap_dims({"nchans":"wavenumber"})
        wavelength=1e+04/ds1.wavenumber
        chkwvn_list=ds1.wavenumber.sel(wavenumber=spectral_range)[ds1.use_flag.sel(wavenumber=spectral_range)==1]
        dates_count+=1
    else:
        print('%s is not existing'%(raddfile0))
        continue


    # Observation lat/lon from exp 0 (baseline)
    rlat0=np.reshape(ds0.Latitude.values,(npts0,nchs0))
    rlon0=np.reshape(ds0.Longitude.values,(npts0,nchs0))
    qcflags0=np.reshape(ds0.QC_Flag.values,(npts0,nchs0))
    #obs0=np.reshape(ds0.Observation.values,(npts0,nchs0))
    sim0=np.reshape(ds0.Simulated_Tb.values,(npts0,nchs0))
    clr0=np.reshape(ds0.Clearsky_Tb.values,(npts0,nchs0))
    varinv0=np.reshape(ds0.Inverse_Observation_Error.values,(npts0,nchs0))
    sim_bc0=np.reshape(ds0.Obs_Minus_Forecast_adjusted.values,(npts0,nchs0))
    sim_nbc0=np.reshape(ds0.Obs_Minus_Forecast_unadjusted.values,(npts0,nchs0))
    obs0=sim_nbc0+sim0
    tmpds0=xa.Dataset({'rlon':(['obsloc'],rlon0[:,0]),
                      'rlat':(['obsloc'],rlat0[:,0]),
                      'qcflag':(['obsloc','wavenumber'],qcflags0),
                       'tb_obs':(['obsloc','wavenumber'],obs0),
                       'tb_sim':(['obsloc','wavenumber'],sim0),
                       'tb_clr':(['obsloc','wavenumber'],clr0),
                      'varinv':(['obsloc','wavenumber'],varinv0),
                      'omb_bc':(['obsloc','wavenumber'],sim_bc0),
                      'omb_nbc':(['obsloc','wavenumber'],sim_nbc0)},
                      coords={'obsloc':np.arange(npts0),
                             'wavenumber':ds0.wavenumber.values})
    tmpds0=tmpds0.sel(wavenumber=chkwvn_list)

    # Observation lat/lon from exp 1 (test)
    rlat1=np.reshape(ds1.Latitude.values,(npts1,nchs1))
    rlon1=np.reshape(ds1.Longitude.values,(npts1,nchs1))
    qcflags1=np.reshape(ds1.QC_Flag.values,(npts1,nchs1))
    #obs1=np.reshape(ds1.Observation.values,(npts1,nchs1))
    sim1=np.reshape(ds1.Simulated_Tb.values,(npts1,nchs1))
    clr1=np.reshape(ds1.Clearsky_Tb.values,(npts1,nchs1))
    varinv1=np.reshape(ds1.Inverse_Observation_Error.values,(npts1,nchs1))
    sim_bc1=np.reshape(ds1.Obs_Minus_Forecast_adjusted.values,(npts1,nchs1))
    sim_nbc1=np.reshape(ds1.Obs_Minus_Forecast_unadjusted.values,(npts1,nchs1))
    obs1=sim_nbc1+sim1
    tmpds1=xa.Dataset({'rlon':(['obsloc'],rlon1[:,0]),
                      'rlat':(['obsloc'],rlat1[:,0]),
                      'qcflag':(['obsloc','wavenumber'],qcflags1),
                       'tb_obs':(['obsloc','wavenumber'],obs1),
                       'tb_sim':(['obsloc','wavenumber'],sim1),
                       'tb_clr':(['obsloc','wavenumber'],clr1),
                      'varinv':(['obsloc','wavenumber'],varinv1),
                      'omb_bc':(['obsloc','wavenumber'],sim_bc1),
                      'omb_nbc':(['obsloc','wavenumber'],sim_nbc1)},
                      coords={'obsloc':np.arange(npts1),
                             'wavenumber':ds1.wavenumber.values})
    tmpds1=tmpds1.sel(wavenumber=chkwvn_list)

    mask0=~np.isnan(tmpds0.rlon)
    mask1=~np.isnan(tmpds1.rlon)

    if (area!='Glb'):
        mask0=(mask0)&((tmpds0.rlon<=maxlon)&(tmpds0.rlon>=minlon)&(tmpds0.rlat<=maxlat)&(tmpds0.rlat>=minlat))
        mask1=(mask1)&((tmpds1.rlon<=maxlon)&(tmpds1.rlon>=minlon)&(tmpds1.rlat<=maxlat)&(tmpds1.rlat>=minlat))

    if (useqc):
        mask0=(mask0)&((tmpds0.qcflag==0)|(tmpds0.qcflag==13))
        mask1=(mask1)&((tmpds1.qcflag==0)|(tmpds1.qcflag==13))

    if (usebc):
        omb0=tmpds0.omb_bc
        omb1=tmpds1.omb_bc
    else:
        omb0=tmpds0.omb_nbc
        omb1=tmpds1.omb_nbc

#    if (date==str(sdate)):
#        omb0_all=omb0
#        ds_all1=tmpds1
#    else:
#        ds_all0=xa.concat((ds_all0,tmpds0),dim='obsloc')
#        ds_all1=xa.concat((ds_all1,tmpds1),dim='obsloc')
#
#total_obscounts=ds_all0.obsloc.size
#ds_all0=ds_all0.assign_coords(obsloc=np.arange(total_obscounts))
#ds_all1=ds_all1.assign_coords(obsloc=np.arange(total_obscounts))
#
#    omf1=xa.where(mask1,np.nan,omf1)
#    omf2=xa.where(mask2,np.nan,omf2)
#    
#    bias[d,0]=np.nanmean(omf1)
#    bias[d,1]=np.nanmean(omf2)
#    rms[d,0]=np.sqrt(np.nanmean(np.square(omf1)))
#    rms[d,1]=np.sqrt(np.nanmean(np.square(omf2)))
#    
#    dates_count+=1
#
#xaxis=np.arange(len(dlist))
#fig,ax=plt.subplots(2,1,sharex=True,figsize=(axe_w,axe_h))
#fig.subplots_adjust(hspace=0.15)
#for a in np.arange(2):
#    ax[a].set_prop_cycle(color=['blue','red'])
#    ax[a].grid(axis='x')
#        
#ax[1].plot_date(dates,bias,'-')
#ax[0].xaxis.set_major_locator(loc)
#ax[0].xaxis.set_major_formatter(formatter)
#ax[0].xaxis.set_tick_params(rotation=30)
#ax[0].set_title('%s %s %.2fµm: AOD>%.2f %s>%.2f'%(
#        sensor,tlstr,wavelength[chkwvlidx],aodmin,aersp,aerfracmin),
#            loc='left')
#ax[0].plot_date(dates,rms,'--')
#ax[0].set_ylabel('RMS [K]')
#ax[1].set_ylabel('BIAS [K]')
#lglist=np.zeros((2,2),dtype='<U30')
#for ex in np.arange(2):
#    lglist[0,ex]=leglist[ex]+'(%8.4f)' %(np.nanmean(rms[:,ex]))
#    lglist[1,ex]=leglist[ex]+'(%8.4f)' %(np.nanmean(bias[:,ex]))
#ax[0].legend(lglist[0,:])
#ax[1].legend(lglist[1,:])
#
#savedir=outpath+'/Innov/1ch'
#if ( not os.path.exists(savedir) ):
#    os.makedirs(savedir)
#
#if (fsave):
#    outname='%s/%s_%s_%s_%s_%s_%s_%s_%s_%.2fµm_BIASRMS.png' %(savedir,area,sensor,tlstr,leglist[0],leglist[1],
#                  qcflg,bcflg,waterflg,wavelength[chkwvlidx])
#    fig.savefig(outname, dpi=quality)
#    plt.close()
    

