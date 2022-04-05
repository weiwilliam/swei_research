# -*- coding: utf-8 -*-
"""
Created on Tue Jun 29 16:06:23 2021

@author: ck102

Aerosol detection based on CADS 3.1 from NWP SAF

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
import numpy as np
import xarray as xa
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mpcrs
import cartopy.crs as ccrs
import setuparea as setarea
from plot_utils import setupax_2dmap, plt_x2y, set_size
from utils import ndate,setup_cmap
from datetime import datetime, timedelta
from matplotlib.dates import (DAILY, DateFormatter,
                              rrulewrapper, RRuleLocator)

tlsize=12 ; lbsize=12
mpl.rc('axes', titlesize=tlsize,labelsize=lbsize)
mpl.rc('xtick',labelsize=lbsize)
mpl.rc('ytick',labelsize=lbsize)
mpl.rc('legend',fontsize='large')
fsave=1 ; ffmt='png' ; ptsize=4
axe_w=6 ; axe_h=3 ; quality=300

# Projection setting
proj=ccrs.PlateCarree(globe=None)

# Plotting setup
sdate=2020061000
edate=2020071018
aertag='Dust'
hint=6
explist=['hazyda_ctrl','hazyda_aero']
expnlist=['CTL','AER']
sensor='iasi_metop-a'
spectral_range=slice(700,1300)
loop='ges' #ges,anl
#if loop=='anl':
#    tlstr='OMA'
#elif loop=='ges':
#    tlstr='OMF'

area='Glb'
minlon, maxlon, minlat, maxlat, crosszero, cyclic=setarea.setarea(area)
print(minlat,maxlat,minlon,maxlon,crosszero,cyclic)
if (area=='Glb'):
   minlon=-180. ; maxlon=180.
cornll=[minlat,maxlat,minlon,maxlon]

cbori='vertical' #vertical, horizontal
if (cbori=='vertical'):
   cb_frac=0.025
   cb_pad=0.06
elif (cbori=='horizontal'):
   cb_frac=0.04
   cb_pad=0.1

# Data path setup
archdir0=rootarch+'/'+explist[0]
archdir1=rootarch+'/'+explist[1]
outpath=rootpath+'/AlbanyWork/Prospectus/Experiments/HazyDA/Images/DiagFiles/rad'

syy=int(str(sdate)[:4]); smm=int(str(sdate)[4:6])
sdd=int(str(sdate)[6:8]); shh=int(str(sdate)[8:10])
eyy=int(str(edate)[:4]); emm=int(str(edate)[4:6])
edd=int(str(edate)[6:8]); ehh=int(str(edate)[8:10])

date1 = datetime(syy,smm,sdd,shh)
date2 = datetime(eyy,emm,edd,ehh)
delta = timedelta(hours=hint)
dates = pd.date_range(start=date1, end=date2, freq=delta)

rule = rrulewrapper(DAILY, byhour=6, interval=4)
loc = RRuleLocator(rule)
formatter = DateFormatter('%Y %h %n %d %Hz')

tnum=0
dlist=[]
cdate=sdate
while (cdate<=edate):
    dlist.append(str(cdate))
    tnum=tnum+1
    cdate=ndate(hint,cdate)

didx=0
for date in dlist:
    didx+=1
    cur_date=dates[didx-1]
    raddfile='diag_'+sensor+'_'+loop+'.'+str(date)+'.nc4'
    infile0=archdir0+'/'+str(date)+'/'+raddfile
    infile1=archdir1+'/'+str(date)+'/'+raddfile

    if (os.path.exists(infile0) and
        os.path.exists(infile1)):
        print('Processing Radfile: %s' %(raddfile))
        ds0=xa.open_dataset(infile0)
        ds1=xa.open_dataset(infile1)
        npts0=int(ds0.nobs.size/ds0.nchans.size)
        npts1=int(ds1.nobs.size/ds1.nchans.size)
        nchs0=ds0.nchans.size
        nchs1=ds1.nchans.size
        ds0=ds0.assign_coords(nuchan=('wavenumber',ds0.wavenumber.data))
        ds0=ds0.swap_dims({"nchans":"wavenumber"})
        ds1=ds1.assign_coords(nuchan=('wavenumber',ds1.wavenumber.data))
        ds1=ds1.swap_dims({"nchans":"wavenumber"}) #replace the dimension of channel by channel indices
        wavelength=1e+04/ds1.wavenumber
        chkwvn_list=ds1.wavenumber.sel(wavenumber=spectral_range)[ds1.use_flag.sel(wavenumber=spectral_range)==1]

        # Observation lat/lon from exp 0 (baseline)
        rlat0=np.reshape(ds0.Latitude.values,(npts0,nchs0))[:,0]
        rlon0=np.reshape(ds0.Longitude.values,(npts0,nchs0))[:,0]
        qcflags0=np.reshape(ds0.QC_Flag.values,(npts0,nchs0))
        # obs0=np.reshape(ds0.Observation.values,(npts0,nchs0))
        sim0=np.reshape(ds0.Simulated_Tb.values,(npts0,nchs0))
        clr0=np.reshape(ds0.Clearsky_Tb.values,(npts0,nchs0))
        varinv0=np.reshape(ds0.Inverse_Observation_Error.values,(npts0,nchs0))
        omb_bc0=np.reshape(ds0.Obs_Minus_Forecast_adjusted.values,(npts0,nchs0))
        omb_nbc0=np.reshape(ds0.Obs_Minus_Forecast_unadjusted.values,(npts0,nchs0))
        obs0=omb_nbc0+sim0
        arrshape=['obsloc','wavenumber']
        tmpds0=xa.Dataset({'rlon':(['obsloc'],rlon0),
                           'rlat':(['obsloc'],rlat0),
                           'qcflag':(arrshape,qcflags0),
                           'tb_obs':(arrshape,obs0),
                           'tb_sim':(arrshape,sim0),
                           'tb_clr':(arrshape,clr0),
                           'varinv':(arrshape,varinv0),
                           'omb_bc':(arrshape,omb_bc0),
                           'omb_nbc':(arrshape,omb_nbc0)},
                           coords={'obsloc':np.arange(npts0),
                                   'wavenumber':ds0.wavenumber.values})
        tmpds0=tmpds0.sel(wavenumber=chkwvn_list)

        # Observation lat/lon from exp 1 (test)
        rlat1=np.reshape(ds1.Latitude.values,(npts1,nchs1))[:,0]
        rlon1=np.reshape(ds1.Longitude.values,(npts1,nchs1))[:,0]
        qcflags1=np.reshape(ds1.QC_Flag.values,(npts1,nchs1))
        # obs1=np.reshape(ds1.Observation.values,(npts1,nchs1))
        sim1=np.reshape(ds1.Simulated_Tb.values,(npts1,nchs1))
        clr1=np.reshape(ds1.Clearsky_Tb.values,(npts1,nchs1))
        varinv1=np.reshape(ds1.Inverse_Observation_Error.values,(npts1,nchs1))
        omb_bc1=np.reshape(ds1.Obs_Minus_Forecast_adjusted.values,(npts1,nchs1))
        omb_nbc1=np.reshape(ds1.Obs_Minus_Forecast_unadjusted.values,(npts1,nchs1))
        obs1=omb_nbc1+sim1
        arrshape=['obsloc','wavenumber']
        tmpds1=xa.Dataset({'rlon':(['obsloc'],rlon1),
                           'rlat':(['obsloc'],rlat1),
                           'qcflag':(arrshape,qcflags1),
                           'tb_obs':(arrshape,obs1),
                           'tb_sim':(arrshape,sim1),
                           'tb_clr':(arrshape,clr1),
                           'varinv':(arrshape,varinv1),
                           'omb_bc':(arrshape,omb_bc1),
                           'omb_nbc':(arrshape,omb_nbc1)},
                           coords={'obsloc':np.arange(npts1),
                                   'wavenumber':ds1.wavenumber.values})
        tmpds1=tmpds1.sel(wavenumber=chkwvn_list)
        
        good_cnts0=np.count_nonzero((tmpds0.qcflag==0),axis=0)
        aero_cnts0=np.count_nonzero((tmpds0.qcflag==13),axis=0)
        good_cnts1=np.count_nonzero((tmpds1.qcflag==0),axis=0)
        aero_cnts1=np.count_nonzero((tmpds1.qcflag==13),axis=0)
        
    else:
        print('%s is not existing'%(raddfile))
        if ('tmpds1' in locals()):
            good_cnts0=xa.zeros_like(tmpds0.tb_obs[0,:])
            aero_cnts0=xa.zeros_like(tmpds0.tb_obs[0,:])
            good_cnts0[:]=np.nan
            aero_cnts0[:]=np.nan
            good_cnts1=xa.zeros_like(tmpds1.tb_obs[0,:])
            aero_cnts1=xa.zeros_like(tmpds1.tb_obs[0,:])
            good_cnts1[:]=np.nan
            aero_cnts1[:]=np.nan
        else:
            print('postpone the starting date')
            continue

    aero_por0=aero_cnts0/(aero_cnts0+good_cnts0)
    aero_por1=aero_cnts1/(aero_cnts1+good_cnts1)
    
    used_cnts=xa.Dataset({'exp0cnts':(['wavenumber'],(good_cnts0+aero_cnts0).data),
                          'exp1cnts':(['wavenumber'],(good_cnts1+aero_cnts1).data),
                          'exp0aerp':(['wavenumber'],(aero_por0).data),
                          'exp1aerp':(['wavenumber'],(aero_por1).data)},
                          coords={'wavenumber':chkwvn_list})
    # Observation lat/lon
    if (date==str(sdate)):
        usedcnts_all=used_cnts
    else:
        usedcnts_all=xa.concat((usedcnts_all,used_cnts),dim='dates')

usedcnts_all=usedcnts_all.assign_coords(dates=('dates',dates))

savedir=outpath+'/'+expnlist[1]+'/datausage/timeseries/'+aertag
if ( not os.path.exists(savedir) ):
    os.makedirs(savedir)
# for chkwvn in [906.25]:
for chkwvn in chkwvn_list:
    qc_chk=usedcnts_all.sel(wavenumber=chkwvn)
    
    fig,ax=plt.subplots()
    set_size(axe_w,axe_h,b=0.2)
    ax.plot_date(dates,qc_chk.exp0cnts,'.-',color='red',linewidth=0.8)
    ax.plot_date(dates,qc_chk.exp1cnts,'.-',color='blue',linewidth=0.8)
    
    ax.xaxis.set_major_locator(loc)
    ax.xaxis.set_major_formatter(formatter)
    ax.xaxis.set_tick_params(rotation=20)
    
    tistr=('%s (%.2f $cm^{-1}$)' %(sensor,chkwvn))
    ax.set_title(tistr,loc='left')
    #ax.set_xlabel('Dates')
    ax.set_ylabel('Number of observations')
    ax.legend(expnlist)
            
    if (fsave):
        fname=('TS_%s_%s_%.2f.%s' %(area,sensor,chkwvn,ffmt))
        fig.savefig(savedir+'/'+fname,dpi=quality)
        plt.close()

savedir=outpath+'/'+expnlist[1]+'/datausage/spec/'+aertag
if ( not os.path.exists(savedir) ):
    os.makedirs(savedir)

wvn=usedcnts_all.wavenumber.data
wvl=1e4/wvn
wvnlb='Wavenumber [$cm^{-1}$]' 
wvllb='Wavelength [µm]'
colorlst=['red','blue']
lstylelst=[' ',' ']
mrklst=['o','^']
lglst=expnlist
cntsyaxlb='Number of observations'
aerpyaxlb='Aerosol-affected [%]'
tistr=''

used0_mean=usedcnts_all.exp0cnts.mean(dim='dates')
used1_mean=usedcnts_all.exp1cnts.mean(dim='dates')
pltda=xa.concat((used0_mean,used1_mean),dim='lines').T

fig,ax=plt.subplots()
set_size(axe_w,axe_h,ax=ax,b=0.25)
plt_x2y(pltda,cntsyaxlb,wvn,wvnlb,wvl,wvllb,colorlst,lstylelst,mrklst,tistr,lglst,0,[],ax=ax)

fname=('Spec_%s_%s_cnts.%s-%s.%s' %(area,sensor,spectral_range.start,spectral_range.stop,ffmt))
if (fsave):
    fig.savefig(savedir+'/'+fname,dpi=quality)
    plt.close()

aerp0_mean=usedcnts_all.exp0aerp.mean(dim='dates')*100.
aerp1_mean=usedcnts_all.exp1aerp.mean(dim='dates')*100.
pltda=xa.concat((aerp0_mean,aerp1_mean),dim='lines').T
        
fig,ax=plt.subplots()
set_size(axe_w,axe_h,ax=ax,b=0.25)
plt_x2y(pltda,aerpyaxlb,wvn,wvnlb,wvl,wvllb,colorlst,lstylelst,mrklst,tistr,lglst,0,[],ax=ax)

fname=('Spec_%s_%s_aerp.%s-%s.%s' %(area,sensor,spectral_range.start,spectral_range.stop,ffmt))
if (fsave):
    fig.savefig(savedir+'/'+fname,dpi=quality)
    plt.close()
