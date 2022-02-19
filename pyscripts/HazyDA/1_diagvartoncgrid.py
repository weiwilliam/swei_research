#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 21:59:35 2019

@author: weiwilliam
"""
import os
import sys
sys.path.append('/data/users/swei/resultcheck/R2O_work/DiagFiles/GSI_DiagPackage/PythonScripts/libs')
import setuparea as setarea
import numpy as np
from datetime import datetime
from datetime import timedelta
import matplotlib.dates as mdates
from matplotlib.dates import (DAILY, DateFormatter,
                              rrulewrapper, RRuleLocator)
import matplotlib.ticker as ticker
import xarray as xa
import pandas as pd
import scipy.stats as ss
import itertools
import multiprocessing as mp
#from multiprocessing import Pool #  Process pool
from multiprocessing import sharedctypes

outpath='/Users/weiwilliam/AlbanyWork/2019/1_R2O/Figures/DiagFiles/Satellite'
fsave=1 ; ffmt='png' 

path='/data/users/swei/archive/nc_DiagFiles'
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
''' ! radiance bias correction terms are as follows:
!  pred(1,:)  = global offset
!  pred(2,:)  = zenith angle predictor, is not used and set to zero now
!  pred(3,:)  = cloud liquid water predictor for clear-sky microwave radiance assimilation
!  pred(4,:)  = square of temperature laps rate predictor
!  pred(5,:)  = temperature laps rate predictor
!  pred(6,:)  = cosinusoidal predictor for SSMI/S ascending/descending bias
!  pred(7,:)  = sinusoidal predictor for SSMI/S
!  pred(8,:)  = emissivity sensitivity predictor for land/sea differences
!  pred(9,:)  = fourth order polynomial of angle bias correction
!  pred(10,:) = third order polynomial of angle bias correction
!  pred(11,:) = second order polynomial of angle bias correction
!  pred(12,:) = first order polynomial of angle bias correction  
!  pred(13,:) = sum of term 9 to 12 
!  pred(14,:) = NSST '''

biasterm=1
biastermname=['Total','GlbOffset','ZenithAng','CLW','LapseRateSq','LapseRate',
              'cos','sin','Emiss','Angle','SST']
deg=2.5

expset=1
if (expset==1):
   explist=np.array(['prctrl','prctrl_anal2'])
   leglist=['CTL','AER']
elif (expset==2):
   explist=np.array(['prctrl','praero'])
   leglist=['CTL','CAER']

sensor='iasi_metop-a'

chkwvl=10.3
colormap=['bwr','bwr','bwr']

sdate=2017080100
edate=2017082818
hint=6

area='Glb'
minlon, maxlon, minlat, maxlat, crosszero=setarea.setarea(area)
print(minlat,maxlat,minlon,maxlon,crosszero)

loop='ges' #ges,anl
if loop=='anl':
    tlstr='OMA'
elif loop=='ges':
    tlstr='OMF'
    
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
delta = timedelta(hours=hint)

dates = pd.date_range(start=date1, end=date2, freq=delta)

xdate2= date2+delta
xdates= mdates.drange(date1, xdate2, delta)

tnum=0
dlist=[]
cdate=sdate
while (cdate<=edate):
    dlist.append(str(cdate))
    tnum=tnum+1
    cdate=ndate(cdate,hint)

d=0
for date in dlist:
    raddfile='diag_'+sensor+'_'+loop+'.'+date+'.nc'
    infile1=path+'/'+explist[0]+'/'+raddfile
    infile2=path+'/'+explist[1]+'/'+raddfile
    
    if (os.path.exists(infile1) and os.path.exists(infile2)):
        print('Processing Radfile: %s' %(raddfile))
        ds1=xa.open_dataset(infile1)
        ds2=xa.open_dataset(infile2)
        if (sensor=='hirs4_n19' or sensor=='hirs4_metop-a' 
            or sensor=='hirs4_metop-b'):
            ds1=ds1.sortby(ds1.wavenumber)
            ds2=ds2.sortby(ds2.wavenumber)
#        print('Load Data Elapsed: %f [s]' %(end-start))
        try:
            chkwvlidx
        except NameError:
            wavelength=1e+04/ds1.wavenumber
            usedchidx=np.where(ds1.iuse_rad==1)[0]
            unusedchidx=np.where(ds1.iuse_rad==-1)[0]
            wvldiff=abs(np.subtract(wavelength[usedchidx],chkwvl))
            chkwvlidx=wvldiff.channel[np.where(wvldiff==wvldiff.min())[0][0]]
            print('Check wavelength: %.2f' %(wavelength[chkwvlidx]))
        #if (usedchidx.size==0):
        #    raise SystemExit('No %s channel be used'%(sensor))
    else:
        print('%s is not existing'%(raddfile))
        d=d+1
        continue

    rlat1=ds1.locinfo[0,:]
    rlon1=ds1.locinfo[1,:]
    rlat2=ds2.locinfo[0,:]
    rlon2=ds2.locinfo[1,:]
    wcov1=ds1.locinfo[10,:]
    wcov2=ds2.locinfo[10,:]
    lcov=ds1.locinfo[11,:]
    if (crosszero):
        rlon1[rlon1>=maxlon]=rlon1[rlon1>=maxlon]-360.
        rlon2[rlon2>=maxlon]=rlon2[rlon2>=maxlon]-360.

    if (biasterm==0):
        pltbias1=ds1.tbcnob[chkwvlidx,:]-ds1.tbc[chkwvlidx,:]
        pltbias2=ds2.tbcnob[chkwvlidx,:]-ds2.tbc[chkwvlidx,:]
    elif (biasterm==9):
        pltbias1=ds1.predbias[chkwvlidx,12,:]
        pltbias2=ds2.predbias[chkwvlidx,12,:]
    elif (biasterm==10):
        pltbias1=ds1.predbias[chkwvlidx,13,:]
        pltbias2=ds2.predbias[chkwvlidx,13,:]
    else:
        pltbias1=ds1.predbias[chkwvlidx,biasterm-1,:]
        pltbias2=ds2.predbias[chkwvlidx,biasterm-1,:]
    #pltbias1=ds1.predbias[chkwvlidx,:,:]
    #pltbias2=ds2.predbias[chkwvlidx,:,:]
    
    try:
        scatda1
    except NameError:
        scatda1=pltbias1 ; lon1=rlon1 ; lat1=rlat1
        scatda2=pltbias2 ; lon2=rlon2 ; lat2=rlat2
    else:
        scatda1=xa.concat((scatda1,pltbias1),dim='obsloc')
        lon1=xa.concat((lon1,rlon1),dim='obsloc')
        lat1=xa.concat((lat1,rlat1),dim='obsloc')
        scatda2=xa.concat((scatda2,pltbias2),dim='obsloc')
        lon2=xa.concat((lon2,rlon2),dim='obsloc')
        lat2=xa.concat((lat2,rlat2),dim='obsloc')
    
    d=d+1

nlon=int(360./deg+1)
nlat=int(180./deg+1)
yblksize=int((nlat-1)/4)
xblksize=int((nlon-1)/8)


grid_lon=np.linspace(0,360,nlon)
grid_lon1=np.zeros((nlon-1),dtype='float')
for x in np.arange(nlon-1):
    grid_lon1[x]=0.5*grid_lon[x]+0.5*grid_lon[x+1]
grid_lat=np.linspace(-90,90,nlat)
grid_lat1=np.zeros((nlat-1),dtype='float')
for y in np.arange(nlat-1):
    grid_lat1[y]=0.5*grid_lat[y]+0.5*grid_lat[y+1]

ave_bias1=np.ctypeslib.as_ctypes(np.zeros((nlat-1,nlon-1),dtype='float'))
sh1=sharedctypes.RawArray(ave_bias1._type_, ave_bias1)

ave_bias2=np.ctypeslib.as_ctypes(np.zeros((nlat-1,nlon-1),dtype='float'))
sh2=sharedctypes.RawArray(ave_bias2._type_, ave_bias2)

var_bias1=np.ctypeslib.as_ctypes(np.zeros((nlat-1,nlon-1),dtype='float'))
sh3=sharedctypes.RawArray(var_bias1._type_, var_bias1)

var_bias2=np.ctypeslib.as_ctypes(np.zeros((nlat-1,nlon-1),dtype='float'))
sh4=sharedctypes.RawArray(var_bias2._type_, ave_bias2)

num_bias1=np.ctypeslib.as_ctypes(np.zeros((nlat-1,nlon-1),dtype='int32'))
sh5=sharedctypes.RawArray(num_bias1._type_, num_bias1)

num_bias2=np.ctypeslib.as_ctypes(np.zeros((nlat-1,nlon-1),dtype='int32'))
sh6=sharedctypes.RawArray(num_bias2._type_, num_bias2)

tstat=np.ctypeslib.as_ctypes(np.zeros((nlat-1,nlon-1),dtype='float'))
sh7=sharedctypes.RawArray(tstat._type_, tstat)

pvalue=np.ctypeslib.as_ctypes(np.zeros((nlat-1,nlon-1),dtype='float'))
sh8=sharedctypes.RawArray(pvalue._type_, pvalue)

def fillthevalue(args): 
    window_y, window_x = args
    tmp1 = np.ctypeslib.as_array(sh1)
    tmp2 = np.ctypeslib.as_array(sh2)
    tmp3 = np.ctypeslib.as_array(sh3)
    tmp4 = np.ctypeslib.as_array(sh4)
    tmp5 = np.ctypeslib.as_array(sh5)
    tmp6 = np.ctypeslib.as_array(sh6)
    tmp7 = np.ctypeslib.as_array(sh7)
    tmp8 = np.ctypeslib.as_array(sh8)

    for y in range(window_y, window_y + yblksize):
        for x in range(window_x, window_x + xblksize):
            mask1=((lon1>=grid_lon[x])&(lon1<=grid_lon[x+1])&
                    (lat1>=grid_lat[y])&(lat1<=grid_lat[y+1]))
            mask2=((lon2>=grid_lon[x])&(lon2<=grid_lon[x+1])&
                    (lat2>=grid_lat[y])&(lat2<=grid_lat[y+1]))
            sctmp1=xa.where(mask1,scatda1,np.nan)
            sctmp2=xa.where(mask2,scatda2,np.nan)
            tmp1[y,x]=np.nanmean(sctmp1)
            tmp2[y,x]=np.nanmean(sctmp2)
            tmp3[y,x]=np.nanvar(sctmp1)
            tmp4[y,x]=np.nanvar(sctmp2)
            tmp5[y,x]=np.count_nonzero(mask1)
            tmp6[y,x]=np.count_nonzero(mask2)
            tmp7[y,x], tmp8[y,x]=ss.ttest_ind(sctmp1[mask1==1],sctmp2[mask2==1])
    print('Finish at'+str(window_y)+' '+str(window_x))

window_idxs = [(i, j) for i, j in
               itertools.product(range(0, nlat-1, yblksize),
                                 range(0, nlon-1, xblksize))]
p = mp.Pool(16)
res = p.map(fillthevalue, window_idxs)
ave_bias1 = np.ctypeslib.as_array(sh1)
ave_bias2 = np.ctypeslib.as_array(sh2)
var_bias1 = np.ctypeslib.as_array(sh3)
var_bias2 = np.ctypeslib.as_array(sh4)
num_bias1 = np.ctypeslib.as_array(sh5)
num_bias2 = np.ctypeslib.as_array(sh6)
tstat = np.ctypeslib.as_array(sh7)
pvalue = np.ctypeslib.as_array(sh8)

ncds=xa.Dataset({'ave_bias1':(['lat','lon'],ave_bias1),
                 'ave_bias2':(['lat','lon'],ave_bias2),
                 'var_bias1':(['lat','lon'],var_bias1),
                 'var_bias2':(['lat','lon'],var_bias2),
                 'num_bias1':(['lat','lon'],num_bias1),
                 'num_bias2':(['lat','lon'],num_bias2),
                 'tstat':(['lat','lon'],tstat),
                 'pvalue':(['lat','lon'],pvalue)},
                coords={'lat':grid_lat1,'lon':grid_lon1},
                attrs={'exp1':leglist[0],'exp2':leglist[1],'start_date':str(sdate),'end_date':str(edate),
                       'cycle':str(hint),'loop':loop,'biasterm':biastermname[biasterm]})
fname='./%s_%s_%.2f_%s_pred%i_%.1fx%.1f_parallel.nc' %(leglist[1],sensor,wavelength[chkwvlidx],loop,biasterm,deg,deg)


ncds.to_netcdf(fname)

p.close()
