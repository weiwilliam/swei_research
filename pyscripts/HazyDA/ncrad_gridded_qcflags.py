#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 21:59:35 2019

@author: weiwilliam
"""
import sys, os, platform
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
import numpy as np
from utils import ndate
from datetime import datetime
from datetime import timedelta
import matplotlib.dates as mdates
import xarray as xa
import pandas as pd
import scipy.stats as ss
import itertools
import multiprocessing as mp
#from multiprocessing import Pool #  Process pool
from multiprocessing import sharedctypes
import warnings

#
warnings.filterwarnings('ignore')

degres=2.5
#degres=1

expset=1
if (expset==1):
   explist=np.array(['hazyda_ctrl','hazyda_aero'])
   leglist=['CTL','AER']
elif (expset==2):
   explist=np.array(['prctrl','praero'])
   leglist=['CTL','CAER']

sensor='iasi_metop-a'

sdate=2020061000
edate=2020071018
hint=6

area='Glb'
minlon, maxlon, minlat, maxlat, crosszero, cyclic=setarea.setarea(area)
print(minlat,maxlat,minlon,maxlon,crosszero,cyclic)

loop='ges' #ges,anl
if loop=='anl':
    tlstr='OMA'
elif loop=='ges':
    tlstr='OMF'

useqc=1
if (useqc):
   qcflg='qc'
else:
   qcflg='noqc'
    
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

outpath=rootpath+'/archive/HazyDA/gridded_diag'
archdir0=rootarch+'/'+explist[0]
archdir1=rootarch+'/'+explist[1]

tnum=0
dlist=[]
cdate=sdate
while (cdate<=edate):
    dlist.append(str(cdate))
    tnum=tnum+1
    cdate=ndate(hint,cdate)

d=0
for date in dlist:
    raddfile='diag_'+sensor+'_'+loop+'.'+date+'.nc4'
    infile0=archdir0+'/'+str(date)+'/'+raddfile
    infile1=archdir1+'/'+str(date)+'/'+raddfile
    
    if (os.path.exists(infile0) and
        os.path.exists(infile1)):
        print('Processing Radfile: %s' %(raddfile),flush=1)
        ds0=xa.open_dataset(infile0)
        ds1=xa.open_dataset(infile1)
        npts0=int(ds0.nobs.size/ds0.nchans.size)
        npts1=int(ds1.nobs.size/ds1.nchans.size)
        nchs0=ds0.nchans.size
        nchs1=ds1.nchans.size
        ds0=ds0.assign_coords(nuchan=('wavenumber',ds0.wavenumber.data))
        ds0=ds0.swap_dims({"nchans":"wavenumber"})
        ds1=ds1.assign_coords(nuchan=('wavenumber',ds1.wavenumber.data))
        ds1=ds1.swap_dims({"nchans":"wavenumber"})
        chkwvn_list=ds0.wavenumber[ds0.use_flag==1]
        if (sensor=='hirs4_n19' or sensor=='hirs4_metop-a' 
            or sensor=='hirs4_metop-b'):
            ds0=ds0.sortby(ds0.wavenumber)
            ds1=ds1.sortby(ds1.wavenumber)
    else:
        print('%s is not existing'%(raddfile))
        d=d+1
        continue

    # Observation lat/lon from exp 0 (test)
    rlat0=np.reshape(ds0.Latitude.values,(npts0,nchs0))[:,0]
    rlon0=np.reshape(ds0.Longitude.values,(npts0,nchs0))[:,0]
    rlon0=(rlon0+180)%360-180
    qcflags0=np.reshape(ds0.QC_Flag.values,(npts0,nchs0))
    omb_bc0=np.reshape(ds0.Obs_Minus_Forecast_adjusted.values,(npts0,nchs0))
    omb_nbc0=np.reshape(ds0.Obs_Minus_Forecast_unadjusted.values,(npts0,nchs0))

    # Observation lat/lon from exp 1 (test)
    rlat1=np.reshape(ds1.Latitude.values,(npts1,nchs1))[:,0]
    rlon1=np.reshape(ds1.Longitude.values,(npts1,nchs1))[:,0]
    rlon1=(rlon1+180)%360-180
    qcflags1=np.reshape(ds1.QC_Flag.values,(npts1,nchs1))
    omb_bc1=np.reshape(ds1.Obs_Minus_Forecast_adjusted.values,(npts1,nchs1))
    omb_nbc1=np.reshape(ds1.Obs_Minus_Forecast_unadjusted.values,(npts1,nchs1))

    if (crosszero):
        rlon1[rlon1>=maxlon]=rlon1[rlon1>=maxlon]-360.
        rlon2[rlon2>=maxlon]=rlon2[rlon2>=maxlon]-360.

    tmpds0=xa.Dataset({'rlon':(['obsloc'],rlon0),
                       'rlat':(['obsloc'],rlat0),
                       'qcflag':(['obsloc','wavenumber'],qcflags0),
                       },
                      coords={'obsloc':np.arange(npts0),
                              'wavenumber':ds0.wavenumber.values})
    tmpds0=tmpds0.sel(wavenumber=chkwvn_list)

    tmpds1=xa.Dataset({'rlon':(['obsloc'],rlon1),
                       'rlat':(['obsloc'],rlat1),
                       'qcflag':(['obsloc','wavenumber'],qcflags1),
                       },
                      coords={'obsloc':np.arange(npts1),
                              'wavenumber':ds0.wavenumber.values})
    tmpds1=tmpds1.sel(wavenumber=chkwvn_list)

    if (date==dlist[0]):
        outds0=tmpds0
        outds1=tmpds1
    else:
        outds0=xa.concat((outds0,tmpds0),dim='obsloc')
        outds1=xa.concat((outds1,tmpds1),dim='obsloc')
    
    d=d+1

qclist=list(set(np.unique(outds0.qcflag))|set(np.unique(outds1.qcflag)))

total_obscounts0=outds0.obsloc.size
outds0=outds0.assign_coords(obsloc=np.arange(total_obscounts0))
total_obscounts1=outds1.obsloc.size
outds1=outds1.assign_coords(obsloc=np.arange(total_obscounts1))
print((total_obscounts0,total_obscounts1),flush=1)

latbin=np.arange(-90,90+0.5*degres,degres)
latgrd=np.arange(-90+0.5*degres,90,degres)
lonbin=np.arange(-180,180+0.5*degres,degres)
longrd=np.arange(-180+0.5*degres,180,degres)

df0=outds0.to_dataframe()
df1=outds1.to_dataframe()

qcname_list=[]
qcidx=0
for qcval in qclist:
    qcname='qc%s'%(int(qcval))
    qcname_list.append(qcname)
    df0_qcfilter=((df0['qcflag']==qcval))
    df1_qcfilter=((df1['qcflag']==qcval))
    outdf0=df0.loc[df0_qcfilter,:]
    outdf1=df1.loc[df1_qcfilter,:]

    outdf0['lat']=pd.cut(outdf0['rlat'],bins=latbin,labels=latgrd)
    outdf0['lon']=pd.cut(outdf0['rlon'],bins=lonbin,labels=longrd)
    outdf1['lat']=pd.cut(outdf1['rlat'],bins=latbin,labels=latgrd)
    outdf1['lon']=pd.cut(outdf1['rlon'],bins=lonbin,labels=longrd)
    grp0 = outdf0.groupby(['wavenumber','lat','lon']).agg({'qcflag':['count']})
    grp1 = outdf1.groupby(['wavenumber','lat','lon']).agg({'qcflag':['count']})

    grp0=grp0.rename(columns={'qcflag':qcname})
    grp1=grp1.rename(columns={'qcflag':qcname})

    if (qcidx==0):
       outgrp0=grp0
       outgrp1=grp1
    else:
       outgrp0=pd.concat((outgrp0,grp0),axis=1)
       outgrp1=pd.concat((outgrp1,grp1),axis=1)
      
    qcidx+=1

grdds0=outgrp0.to_xarray()
grdds1=outgrp1.to_xarray()

for var in qcname_list:
   for stats in ['count']:
      newname='%s_%s'%(var,stats)
      grdds0=grdds0.rename({(var,stats):(newname)})
      grdds1=grdds1.rename({(var,stats):(newname)})
#
fname0='%s/%s_%s_%s_qcflags_%.1fx%.1f.%s_%s.nc' %(outpath,leglist[0],sensor,loop,degres,degres,sdate,edate)
print(fname0,flush=1)
grdds0.to_netcdf(fname0)

fname1='%s/%s_%s_%s_qcflags_%.1fx%.1f.%s_%s.nc' %(outpath,leglist[1],sensor,loop,degres,degres,sdate,edate)
print(fname1,flush=1)
grdds1.to_netcdf(fname1)
