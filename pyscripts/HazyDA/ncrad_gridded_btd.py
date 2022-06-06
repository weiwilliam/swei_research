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

sensorlist=['airs_aqua','amsua_aqua','amsua_metop-a','amsua_n15','amsua_n18',
            'amsua_n19','atms_npp','avhrr_metop-a','avhrr_n18','cris_npp','gmi_gpm',
            'hirs4_metop-a','hirs4_metop-b','hirs4_n19','iasi_metop-a','iasi_metop-b',
            'mhs_metop-a','mhs_metop-b','mhs_n18','mhs_n19','saphir_meghat',
            'seviri_m08','seviri_m10','sndrd1_g15','sndrd2_g15','sndrd3_g15',
            'sndrd4_g15','ssmis_f17','ssmis_f18']
#hsensorlist=['airs_aqua','iasi_metop-a','iasi_metop-b','cris_npp']
#lsensor1list=['hirs4_metop-a','hirs4_metop-b','hirs4_n19']
#lsensor2list=['sndrd1_g15','sndrd2_g15','sndrd3_g15','sndrd4_g15']
#lsensor3list=['avhrr_metop-a','avhrr_n18','seviri_m08','seviri_m10']

degres=2.5
#degres=1

explist=['hazyda_aero_sea']
leglist=['AERS']

sensor='iasi_metop-a'

sdate=2020060106
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

latbin=np.arange(-90,90+0.5*degres,degres)
latgrd=np.arange(-90+0.5*degres,90,degres)
lonbin=np.arange(-180,180+0.5*degres,degres)
longrd=np.arange(-180+0.5*degres,180,degres)

tnum=0
dlist=[]
cdate=sdate
while (cdate<=edate):
    dlist.append(str(cdate))
    tnum=tnum+1
    cdate=ndate(hint,cdate)

idx=0
for date in dlist:
    raddfile='diag_'+sensor+'_'+loop+'.'+date+'.nc4'
    infile0=archdir0+'/'+str(date)+'/'+raddfile
    
    if (os.path.exists(infile0)):
        print('Processing Radfile: %s' %(raddfile),flush=1)
        ds0=xa.open_dataset(infile0)
        npts0=int(ds0.nobs.size/ds0.nchans.size)
        nchs0=ds0.nchans.size
        chkwvn_list=ds0.wavenumber[ds0.use_flag==1]
        if (sensor=='hirs4_n19' or sensor=='hirs4_metop-a' 
            or sensor=='hirs4_metop-b'):
            ds0=ds0.sortby(ds0.wavenumber)
            ds1=ds1.sortby(ds1.wavenumber)
        if (idx==0):
           avaldates=dates[idx]
        else:
           avaldates=np.append(avaldates,dates[idx])
        idx+=1
    else:
        print('%s is not existing'%(raddfile))
        idx+=1
        continue

    # Observation lat/lon from exp 0 (test)
    rlat0=np.reshape(ds0.Latitude.values,(npts0,nchs0))[:,0]
    rlon0=np.reshape(ds0.Longitude.values,(npts0,nchs0))[:,0]
    rlon0=(rlon0+180)%360-180
    qcflags0=np.reshape(ds0.QC_Flag.values,(npts0,nchs0))
    sim0=np.reshape(ds0.Simulated_Tb.values,(npts0,nchs0))
    clr0=np.reshape(ds0.Clearsky_Tb.values,(npts0,nchs0))
    btd0=sim0-clr0

    tmpds0=xa.Dataset({'rlon':(['obsloc'],rlon0),
                       'rlat':(['obsloc'],rlat0),
                       'qcflag':(['obsloc','wavenumber'],qcflags0),
                       'btd':(['obsloc','wavenumber'],btd0),
                       },
                      coords={'obsloc':np.arange(npts0),
                              'wavenumber':ds0.wavenumber.values})
    tmpds0=tmpds0.sel(wavenumber=chkwvn_list)

    tmpdf0=tmpds0.to_dataframe()
    if (useqc):
       tmpdf0_qcfilter=((tmpdf0['qcflag']==0.0)|(tmpdf0['qcflag']==13.0))
       tmpoutdf0=tmpdf0.loc[tmpdf0_qcfilter,:]
    else:
       tmpoutdf0=df0
    tmpoutdf0=tmpoutdf0.reset_index()
    tmpoutdf0['lat']=pd.cut(tmpoutdf0['rlat'],bins=latbin,labels=latgrd)
    tmpoutdf0['lon']=pd.cut(tmpoutdf0['rlon'],bins=lonbin,labels=longrd)
    
    tmpgrp0  =tmpoutdf0.groupby(['wavenumber','lat','lon']).agg({'btd':['mean','count','var','max','min']})
    tmpgrdds0=tmpgrp0.to_xarray()
    
    for stats in ['mean','count','var','max','min']:
       newname='btd_%s'%(stats)
       tmpgrdds0=tmpgrdds0.rename({('btd',stats):(newname)})

    if (date==dlist[0]):
        outds0=tmpds0
        tsgrd0=tmpgrdds0
    else:
        outds0=xa.concat((outds0,tmpds0),dim='obsloc')
        tsgrd0=xa.concat((tsgrd0,tmpgrdds0),dim='time')

tsgrd0=tsgrd0.assign_coords({'time':avaldates})

fname0='%s/%s_%s_%s_%s_btd_%.1fx%.1f.time.%s_%s.nc' %(outpath,leglist[0],sensor,loop,qcflg,degres,degres,sdate,edate)
print(fname0,flush=1)
tsgrd0.to_netcdf(fname0)

print('Processing gridded data for whole period',flush=1)

total_obscounts0=outds0.obsloc.size
outds0=outds0.assign_coords(obsloc=np.arange(total_obscounts0))

df0=outds0.to_dataframe()
if (useqc):
   df0_qcfilter=((df0['qcflag']==0.0))
   outdf0=df0.loc[df0_qcfilter,:]
else:
   outdf0=df0
outdf0['lat']=pd.cut(outdf0['rlat'],bins=latbin,labels=latgrd)
outdf0['lon']=pd.cut(outdf0['rlon'],bins=lonbin,labels=longrd)

grp0 = outdf0.groupby(['wavenumber','lat','lon']).agg({'btd':['mean','count','var','max','min']})
grdds0=grp0.to_xarray()

for stats in ['mean','count','var','max','min']:
   newname='btd_%s'%(stats)
   grdds0=grdds0.rename({('btd',stats):(newname)})

fname1='%s/%s_%s_%s_%s_btd_%.1fx%.1f.mean.%s_%s.nc' %(outpath,leglist[0],sensor,loop,qcflg,degres,degres,sdate,edate)
print(fname1,flush=1)
grdds0.to_netcdf(fname1)
