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
    exthdarch='D:\ResearchData'
    rootgit='F:\GitHub\swei_research'
sys.path.append(rootgit+'/pyscripts/functions')
import setuparea as setarea
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
import time
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

outpath=rootpath+'/AlbanyWork/R2O_Project/images/DiagFiles/Satellite'
path=exthdarch+'/2_R2O/nc_DiagFiles'

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

explist=np.array(['prctrl','praero'])
leglist=['CTL','AER']
sensor='iasi_metop-a'

chkwvl=13.
colormap=['bwr','bwr','bwr']

sdate=2017080100
edate=2017082818
hint=6

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
useqc=0
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
    
useaod=0
aersp='dust' #: dust, seas, sulf, carb
if (useaod):
    aodmin=0.3
    aerfracmin=0.65
else:
    aodmin=0.
    aerfracmin=0.

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

rule = rrulewrapper(DAILY, byhour=6, interval=5)
loc = RRuleLocator(rule)
formatter = DateFormatter('%Y%h %n %d %Hz')

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
    infile1=path+'/'+explist[0]+'/'+date+'/'+raddfile
    infile2=path+'/'+explist[1]+'/'+date+'/'+raddfile
    
    if (os.path.exists(infile1) and os.path.exists(infile2)):
        print('Processing Radfile: %s' %(raddfile))
        start=time.time()
        ds1=xa.open_dataset(infile1)
        ds2=xa.open_dataset(infile2)
        if (sensor=='hirs4_n19' or sensor=='hirs4_metop-a' 
            or sensor=='hirs4_metop-b'):
            ds1=ds1.sortby(ds1.wavenumber)
            ds2=ds2.sortby(ds2.wavenumber)
        end=time.time()
#        print('Load Data Elapsed: %f [s]' %(end-start))
        try:
            chkwvlidx
        except NameError:
            wavelength=1e+04/ds1.wavenumber
            usedchidx=np.where(ds1.iuse_rad==1)[0]
            unusedchidx=np.where(ds1.iuse_rad==-1)[0]
            wvldiff=abs(np.subtract(wavelength[usedchidx],chkwvl))
            chkwvlidx=wvldiff.channel[np.where(wvldiff==wvldiff.min())[0][0]]
            bias=np.zeros((tnum,2),dtype='float'); bias[:,:]=np.nan
            rms=np.zeros((tnum,2),dtype='float'); rms[:,:]=np.nan
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

    start=time.time()
    if (area != 'Glb'):
        areaflg1=np.zeros((ds1.obsloc.size),dtype='int')
        areaflg1[rlat1<=maxlat]=1 ; areaflg1[rlat1<=minlat]=0
        areaflg1[rlon1>=maxlon]=0 ; areaflg1[rlon1<=minlon]=0
        areaidx1=np.where(areaflg1==1)[0]
        areaflg2=np.zeros((ds2.obsloc.size),dtype='int')
        areaflg2[rlat2<=maxlat]=1 ; areaflg2[rlat2<=minlat]=0
        areaflg2[rlon2>=maxlon]=0 ; areaflg2[rlon2<=minlon]=0
        areaidx2=np.where(areaflg2==1)[0]
    else:
        areaidx1=ds1.obsloc
        areaidx2=ds2.obsloc
    end=time.time()
#    print('Get %s Areaindex Elapsed: %f [s]' %(area,end-start))
    
    if (aersp=='dust'):
        aerfrac1=ds1.dufrac
        aerfrac2=ds2.dufrac
    if (aersp=='seas'):
        aerfrac1=ds1.ssfrac
        aerfrac2=ds2.ssfrac
    if (aersp=='sulf'):
        aerfrac1=ds1.sufrac
        aerfrac2=ds2.sufrac
    if (aersp=='carbon'):
        aerfrac1=ds1.ocfrac+ds1.bcfrac
        aerfrac2=ds2.ocfrac+ds2.bcfrac
    
    mask1=(ds1.aod<aodmin)|(aerfrac1<aerfracmin)
    mask2=(ds2.aod<aodmin)|(aerfrac2<aerfracmin)
            
    if (area!='Glb'):
        mask1=(mask1)|(areaflg1==0)
        mask2=(mask2)|(areaflg2==0)
        
    if (useqc):
        mask1=(mask1)|(abs(ds1.qcflag[chkwvlidx,:])!=0)
        mask2=(mask2)|(abs(ds2.qcflag[chkwvlidx,:])!=0)
    
    if (wateronly):
        mask1=(mask1)|(wcov1!=1.)
        mask2=(mask2)|(wcov2!=1.)
    
    if (usebc):
        omf1=ds1.tbc[chkwvlidx,:]
        omf2=ds2.tbc[chkwvlidx,:]
    else:
        omf1=ds1.tbcnob[chkwvlidx,:]
        omf2=ds2.tbcnob[chkwvlidx,:]
        
    omf1=xa.where(mask1,np.nan,omf1)
    omf2=xa.where(mask2,np.nan,omf2)
    
    bias[d,0]=np.nanmean(omf1)
    bias[d,1]=np.nanmean(omf2)
    rms[d,0]=np.sqrt(np.nanmean(np.square(omf1)))
    rms[d,1]=np.sqrt(np.nanmean(np.square(omf2)))
    
    d=d+1

xaxis=np.arange(len(dlist))
fig,ax=plt.subplots(2,1,sharex=True,figsize=(axe_w,axe_h))
fig.subplots_adjust(hspace=0.15)
for a in np.arange(2):
    ax[a].set_prop_cycle(color=['blue','red'])
    ax[a].grid(axis='x')
        
ax[1].plot_date(dates,bias,'-')
ax[0].xaxis.set_major_locator(loc)
ax[0].xaxis.set_major_formatter(formatter)
ax[0].xaxis.set_tick_params(rotation=30)
ax[0].set_title('%s %s %.2fµm: AOD>%.2f %s>%.2f'%(
        sensor,tlstr,wavelength[chkwvlidx],aodmin,aersp,aerfracmin),
            loc='left')
ax[0].plot_date(dates,rms,'--')
ax[0].set_ylabel('RMS [K]')
ax[1].set_ylabel('BIAS [K]')
lglist=np.zeros((2,2),dtype='<U30')
for ex in np.arange(2):
    lglist[0,ex]=leglist[ex]+'(%8.4f)' %(np.nanmean(rms[:,ex]))
    lglist[1,ex]=leglist[ex]+'(%8.4f)' %(np.nanmean(bias[:,ex]))
ax[0].legend(lglist[0,:])
ax[1].legend(lglist[1,:])

savedir=outpath+'/Innov/1ch'
if ( not os.path.exists(savedir) ):
    os.makedirs(savedir)

if (fsave):
    outname='%s/%s_%s_%s_%s_%s_%s_%s_%s_%.2fµm_BIASRMS.png' %(savedir,area,sensor,tlstr,leglist[0],leglist[1],
                  qcflg,bcflg,waterflg,wavelength[chkwvlidx])
    fig.savefig(outname, dpi=quality)
    plt.close()
    

