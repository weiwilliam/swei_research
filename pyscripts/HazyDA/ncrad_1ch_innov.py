#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 21:59:35 2019

@author: weiwilliam
"""
import os, sys, platform
machine='S4'
if (machine=='MBP'):
    rootpath='/Users/weiwilliam'
    rootarch='/Volumes/WD2TB/ResearchData'
elif (machine=='Desktop'):
    rootpath='F:\GoogleDrive_NCU\Albany'
    rootarch='F:\ResearchData'
    rootgit='F:\GitHub\swei_research'
elif (machine=='S4'):
    rootpath='/data/users/swei'
    rootarch='/scratch/users/swei/ncdiag'
    rootgit='/home/swei/research'
elif (machine=='Hera'):
    rootpath='/scratch2/BMC/gsd-fv3-dev/Shih-wei.Wei'
    rootarch='/scratch2/BMC/gsd-fv3-dev/Shih-wei.Wei/ResearchData'
    rootgit='/home/Shih-wei.Wei/research'
elif (machine=='Cheyenne'):
    rootpath='/glade/work/swei/output/images'
    rootarch='/scratch2/BMC/gsd-fv3-dev/Shih-wei.Wei/ResearchData'
    rootgit='/glade/u/home/swei/research'
sys.path.append(rootgit+'/pyscripts/functions')
import setuparea as setarea
from utils import ndate 
from plot_utils import set_size
from gsi_ncdiag import read_rad_ncdiag
from datetime import datetime
from datetime import timedelta
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.dates import (DAILY, DateFormatter,
                              rrulewrapper, RRuleLocator)
import matplotlib.ticker as ticker
import numpy as np
import xarray as xa
import pandas as pd

tlsize=10 ; lbsize=8
mpl.rc('axes', titlesize=tlsize,labelsize=lbsize)
mpl.rc('xtick',labelsize=lbsize)
mpl.rc('ytick',labelsize=lbsize)
mpl.rc('legend',fontsize='small')
fsave=1 ; ffmt='png' ; ptsize=4
axe_w=6 ; axe_h=2.7 ; quality=300
diagsuffix='nc4'

outputpath=rootpath+'/AlbanyWork/Prospectus/Experiments/HazyDA/Images'
inputpath=rootarch

explist=np.array(['hazyda_ctrl','hazyda_aero'])
expnlist=['CTL','AER']
sensor='iasi_metop-a'
units='K'

sdate=2020061000
edate=2020071018
hint=6
chkwvn=962.5

area='Glb'
minlon, maxlon, minlat, maxlat, crosszero, cyclic=setarea.setarea(area)
print(area,minlat,maxlat,minlon,maxlon,crosszero,cyclic)

loop='ges' #ges,anl
usebc=1
# -2: no qc, -1: qc, other: flag
useqc=-1
# -1: exp self msk, 0: use first exp msk, 1: use second exp msk 2: intersection of both exp (useqc=0)
usemsk=3

if loop=='anl':
   lpstr='OMA'
elif loop=='ges':
   lpstr='OMF'

if (usebc):
    bcflg='bc'
    bc_ylbstr='w/ BC'
else:
    bcflg='nobc'
    bc_ylbstr='w/o BC'

if (useqc==-2):
    qcflg='noqc'
elif (useqc==-1):
    qcflg='qc'
else:
    qcflg='qc%s'%(useqc)

rmflier=1
wateronly=0
if (wateronly):
    waterflg='water'
else:
    waterflg='all'
if (usemsk==-1):
    mskflg='omsk'
elif (usemsk==0):
    mskflg='msk0'
elif (usemsk==1):
    mskflg='msk1'
elif (usemsk==2):
    mskflg='imsk'
elif (usemsk==3):
    mskflg='iaer'

imgsavpath=outputpath+'/DiagFiles/rad/1ch/omb/'+area
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

idx=0
for date in dlist:
    raddfile='diag_'+sensor+'_'+loop+'.'+date+'.'+diagsuffix
    infile0=inputpath+'/'+explist[0]+'/'+date+'/'+raddfile
    infile1=inputpath+'/'+explist[1]+'/'+date+'/'+raddfile
    
    if (os.path.exists(infile0) and
        os.path.exists(infile1) ):
        print('Processing Radfile: %s' %(raddfile),flush=1)
        ds0=xa.open_dataset(infile0)
        npts0=int(ds0.nobs.size/ds0.nchans.size)
        nchs0=ds0.nchans.size
        if (idx==0):
           avaldates=dates[idx]
        else:
           avaldates=np.append(avaldates,dates[idx])
        idx+=1
    else:
        print('%s is not existing'%(raddfile),flush=1)
        idx+=1
        continue

    tmpds0=read_rad_ncdiag(infile0,chkwvn=chkwvn)
    tmpds1=read_rad_ncdiag(infile1,chkwvn=chkwvn)
#
    if (date==str(sdate)):
        ds_all0=tmpds0
        ds_all1=tmpds1
    else:
        ds_all0=xa.concat((ds_all0,tmpds0),dim='time')
        ds_all1=xa.concat((ds_all1,tmpds1),dim='time')
ds_all0.assign_coords({'time':avaldates})
ds_all1.assign_coords({'time':avaldates})
#

total_obscounts=ds_all0.obsloc.size
ds_all0=ds_all0.assign_coords(obsloc=np.arange(total_obscounts))
ds_all1=ds_all1.assign_coords(obsloc=np.arange(total_obscounts))
#
mask0=~np.isnan(ds_all0.rlon)
mask1=~np.isnan(ds_all1.rlon)

if (area!='Glb'):
    mask0=(mask0)&((ds_all0.rlon<=maxlon)&(ds_all0.rlon>=minlon)&(ds_all0.rlat<=maxlat)&(ds_all0.rlat>=minlat))
    mask1=(mask1)&((ds_all1.rlon<=maxlon)&(ds_all1.rlon>=minlon)&(ds_all1.rlat<=maxlat)&(ds_all1.rlat>=minlat))

if (useqc==-2):
    pass
elif (useqc==-1):
    mask0=(mask0)&((ds_all0.qcflag==0)|(ds_all0.qcflag==13))
    mask1=(mask1)&((ds_all1.qcflag==0)|(ds_all1.qcflag==13))
else:
    mask0=(mask0)&((ds_all0.qcflag==useqc))
    mask1=(mask1)&((ds_all1.qcflag==useqc))
#
if (usebc):
    omb0=ds_all0.omb_bc
    omb1=ds_all1.omb_bc
else:
    omb0=ds_all0.omb_nbc
    omb1=ds_all1.omb_nbc

if (usemsk==-1):
    pltmsk0=mask0
    pltmsk1=mask1
elif (usemsk==0):
    pltmsk0=mask0
    pltmsk1=mask0
elif (usemsk==1):
    pltmsk0=mask1
    pltmsk1=mask1
elif (usemsk==2):
    pltmsk0=mask0&mask1
    pltmsk1=mask0&mask1
elif (usemsk==3):
    pltmsk0=mask0&mask1&(ds_all1.qcflag==13)
    pltmsk1=mask0&mask1&(ds_all1.qcflag==13)

print('plot counts: %s, %s' %(np.count_nonzero(pltmsk0),np.count_nonzero(pltmsk1)))

omb0=xa.where(pltmsk0,omb0,np.nan)
omb1=xa.where(pltmsk1,omb1,np.nan)

omb_mean0=omb0.mean(dim='obsloc',skipna=1)
omb_mean1=omb1.mean(dim='obsloc',skipna=1)

omb_rms0=np.sqrt((omb0*omb0).mean(dim='obsloc',skipna=1))
omb_rms1=np.sqrt((omb1*omb1).mean(dim='obsloc',skipna=1))

omb=xa.concat((omb_mean0,omb_mean1),dim='exp')
rms=xa.concat((omb_rms0,omb_rms1),dim='exp')
pltmean=omb.data.swapaxes(0,1)
pltrms=rms.data.swapaxes(0,1)

fig,ax=plt.subplots(2,1,sharex=True,figsize=(9,3.8))
fig.subplots_adjust(hspace=0.1)
for a in np.arange(2):
    ax[a].set_prop_cycle(color=['blue','red'])
    ax[a].grid(axis='x')

ax[1].plot_date(avaldates,pltmean,'-o',ms=2)
ax[0].xaxis.set_major_locator(loc)
ax[0].xaxis.set_major_formatter(formatter)
ax[0].xaxis.set_tick_params(rotation=30, labelsize=10)
ax[0].plot_date(avaldates,pltrms,'-o',ms=2)
ax[0].set_ylabel('RMS %s [%s]' %(lpstr,units))
ax[1].set_ylabel('Mean %s [%s]'%(lpstr,units))
lglist=np.zeros((2,2),dtype='<U30')
for ex in np.arange(2):
    lglist[0,ex]=expnlist[ex]+'(%.2f)' %(np.nanmean(pltrms[:,ex]))
    lglist[1,ex]=expnlist[ex]+'(%.2f)' %(np.nanmean(pltmean[:,ex]))
ax[0].legend(lglist[0,:])
ax[1].legend(lglist[1,:])
tistr='%s %.2f $\mathrm{cm^{-1}}$'%(sensor,chkwvn)
ax[0].set_title(tistr,loc='left')
#
if (fsave):
    outname=('%s/%s_%s_%s_%s_%s_%s_%s_%s_%s_%.2f_BIASRMS.%s_%s.png'
             %(imgsavpath,area,sensor,lpstr,
               expnlist[0],expnlist[1],bcflg,qcflg,mskflg,waterflg,chkwvn,sdate,edate))
    print(outname,flush=1)
    fig.savefig(outname, dpi=quality)
    plt.close()
    

