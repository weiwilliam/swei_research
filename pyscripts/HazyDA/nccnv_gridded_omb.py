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
import warnings

#
warnings.filterwarnings('ignore')

degres=2.5
#degres=1
statslist=['mean','count','var']#,'max','min']

explist=['hazyda_aero']
leglist=['AER']

cnvvar='q'
if (cnvvar in ['u','v']):
   cnvtype='uv'
else:
   cnvtype=cnvvar

if (cnvtype in ['ps','sst','tcp']):
    is_sfc=1
elif (cnvtype in ['gps','q','t','uv']):
    is_sfc=0

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

if (is_sfc):
   typeflg='sfc'
else:
   typeflg='upa'
    
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

ptop=np.array((1000.,900.,800.,600.,400.,300.,250.,200.,150.,100.,50.,0.))
pbot=np.array((1200.,1000.,900.,800.,600.,400.,300.,250.,200.,150.,100.,50.))
levgrd=np.sort(0.5*(ptop+pbot))
levbin=np.union1d(ptop,pbot)

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
    #print('%s Start' %(datetime.now().strftime('%c')),flush=1)
    cnvdfile='diag_conv_'+cnvtype+'_'+loop+'.'+date+'.nc4'
    infile0=archdir0+'/'+str(date)+'/'+cnvdfile
    
    if (os.path.exists(infile0)):
        print('Processing Cnvfile: %s' %(cnvdfile),flush=1)
        ds0=xa.open_dataset(infile0)
        npts0=ds0.nobs.data
        if (idx==0):
           avaldates=dates[idx]
        else:
           avaldates=np.append(avaldates,dates[idx])
        idx+=1
    else:
        print('%s is not existing'%(cnvdfile))
        idx+=1
        continue

    # Observation lat/lon from exp 0 (test)
    rlat0=ds0.Latitude.data
    rlon0=ds0.Longitude.data
    rlon0=(rlon0+180)%360-180
    pres0=ds0.Pressure.data
    qcflags0=ds0.Analysis_Use_Flag.data
    obstype0=ds0.Observation_Type.data
    if (cnvtype=='uv'):
       if (cnvvar=='u'):
          omb_bc0 =ds0.u_Obs_Minus_Forecast_adjusted.data
          omb_nbc0=ds0.u_Obs_Minus_Forecast_unadjusted.data
       elif (cnvvar=='v'):
          omb_bc0 =ds0.v_Obs_Minus_Forecast_adjusted.data
          omb_nbc0=ds0.v_Obs_Minus_Forecast_unadjusted.data
    else:
       omb_bc0 =ds0.Obs_Minus_Forecast_adjusted.data
       omb_nbc0=ds0.Obs_Minus_Forecast_unadjusted.data
    varlist=['omb_bc','omb_nbc']

    tmpds0=xa.Dataset({'rlon':(['obsloc'],rlon0),
                       'rlat':(['obsloc'],rlat0),
                       'pres':(['obsloc'],pres0),
                       'qcflag':(['obsloc'],qcflags0),
                       'obstype':(['obsloc'],obstype0),
                       'omb_bc':(['obsloc'],omb_bc0),
                       'omb_nbc':(['obsloc'],omb_nbc0),
                       },
                      coords={'obsloc':npts0})

    tmpdf0=tmpds0.to_dataframe()
    if (useqc):
       tmpdf0_qcfilter=(tmpdf0['qcflag']==1.)
       tmpoutdf0=tmpdf0.loc[tmpdf0_qcfilter,:]
    else:
       tmpoutdf0=tmpdf0
    tmpoutdf0=tmpoutdf0.reset_index()
    tmpoutdf0['lat']=pd.cut(tmpoutdf0['rlat'],bins=latbin,labels=latgrd)
    tmpoutdf0['lon']=pd.cut(tmpoutdf0['rlon'],bins=lonbin,labels=longrd)
    if (is_sfc):
       tmpgrp0 = tmpoutdf0.groupby(['obstype','lat','lon']).agg({'omb_bc':statslist,
                                                                 'omb_nbc':statslist})
    else:
       tmpoutdf0['lev']=pd.cut(tmpoutdf0['pres'],bins=levbin,labels=levgrd)
       tmpgrp0 = tmpoutdf0.groupby(['obstype','lev','lat','lon']).agg({'omb_bc':statslist,
                                                                       'omb_nbc':statslist})

    tmpgrdds0=tmpgrp0.to_xarray()

    for var in varlist:
       for stats in statslist:
          newname='%s_%s'%(var,stats)
          tmpgrdds0=tmpgrdds0.rename({(var,stats):(newname)})

    if (date==dlist[0]):
        outds0=tmpds0
        tsgrd0=tmpgrdds0
    else:
        outds0=xa.concat((outds0,tmpds0),dim='obsloc')
        tsgrd0=xa.concat((tsgrd0,tmpgrdds0),dim='time')
print('%s End' %(datetime.now().strftime('%c')),flush=1)

tsgrd0=tsgrd0.assign_coords({'time':avaldates})

fname0='%s/%s_%s_%s_%s_omb_%.1fx%.1f.time.%s_%s.nc' %(outpath,leglist[0],cnvvar,loop,qcflg,degres,degres,sdate,edate)
print(fname0,flush=1)
tsgrd0.to_netcdf(fname0)

print('%s Processing gridded data for whole period'%(datetime.now().strftime('%c')),flush=1)

total_obscounts0=outds0.obsloc.size
outds0=outds0.assign_coords(obsloc=np.arange(total_obscounts0))

if (is_sfc):
   base_vars=['rlon','rlat','qcflag','obstype']
else:
   base_vars=['rlon','rlat','pres','qcflag','obstype']

#vidx=0
#tmp_vars=base_vars[:]
#tmp_vars.append(var)
#outds0=outds0[tmp_vars]

df0=outds0.to_dataframe()
if (useqc):
   df0_qcfilter=(df0['qcflag']==1.)
   outdf0=df0.loc[df0_qcfilter,:]
else:
   outdf0=df0
outdf0['lat']=pd.cut(outdf0['rlat'],bins=latbin,labels=latgrd)
outdf0['lon']=pd.cut(outdf0['rlon'],bins=lonbin,labels=longrd)

if (is_sfc):
   grp0 = outdf0.groupby(['obstype','lat','lon']).agg({'omb_bc':statslist,
                                                       'omb_nbc':statslist})
else:
   outdf0['lev']=pd.cut(outdf0['pres'],bins=levbin,labels=levgrd)
   grp0 = outdf0.groupby(['obstype','lev','lat','lon']).agg({'omb_bc':statslist,
                                                             'omb_nbc':statslist})

grdds0=grp0.to_xarray()

for var in varlist:
   for stats in statslist:
      newname='%s_%s'%(var,stats)
      grdds0=grdds0.rename({(var,stats):(newname)})
#   if (vidx==0):
#      grdds0=tmpgrdds0
#   else:
#      grdds0=xa.merge((grdds0,tmpgrdds0))

#   vidx+=1

print('%s End Processing whole period'%(datetime.now().strftime('%c')),flush=1)

fname1='%s/%s_%s_%s_%s_omb_%.1fx%.1f.mean.%s_%s.nc' %(outpath,leglist[0],cnvvar,loop,qcflg,degres,degres,sdate,edate)
print(fname1,flush=1)
grdds0.to_netcdf(fname1)
