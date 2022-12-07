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
    rootarch='/data/users/swei/archive/nc_DiagFiles'
    rootgit='/home/swei/research'
elif (machine=='Hera'):
    rootpath='/scratch2/BMC/gsd-fv3-dev/Shih-wei.Wei'
    rootarch='/scratch2/BMC/gsd-fv3-dev/Shih-wei.Wei/ResearchData'
    rootgit='/home/Shih-wei.Wei/research'
elif (machine=='Cheyenne'):
    rootpath='/glade/work/swei/output/images'
    rootarch='/scratch2/BMC/gsd-fv3-dev/Shih-wei.Wei/ResearchData'
    rootgit='/glade/u/home/swei/research'
import numpy as np
import xarray as xa
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mpcrs
import cartopy.crs as ccrs
import setuparea as setarea
from plot_utils import setupax_2dmap, plt_2exps_x2y, set_size
from utils import ndate,setup_cmap
from datetime import datetime, timedelta
from matplotlib.dates import (DAILY, DateFormatter,
                              rrulewrapper, RRuleLocator)
from gsi_ncdiag import read_rad_ncdiag

tlsize=12 ; lbsize=12
mpl.rc('axes', titlesize=tlsize,labelsize=lbsize)
mpl.rc('xtick',labelsize=lbsize)
mpl.rc('ytick',labelsize=lbsize)
mpl.rc('legend',fontsize='large')
fsave=1 ; ffmt='png' ; ptsize=4
axe_w=6 ; axe_h=3 ; quality=300
minussign=u'\u2212'

# Projection setting
proj=ccrs.PlateCarree(globe=None)

# Plotting setup
sdate=2020061000
edate=2020071018
hint=6
explist=['hazyda_ctrl','hazyda_aero']
expnlist=['CTL','AER']
sensor='iasi_metop-a'
spectral_range=slice(750,1200)
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

ts_savedir=outpath+'/datausage/ts/'+area
if ( not os.path.exists(ts_savedir) ):
    os.makedirs(ts_savedir)
sp_savedir=outpath+'/datausage/sp/'+area
if ( not os.path.exists(sp_savedir) ):
    os.makedirs(sp_savedir)

date1 = pd.to_datetime(sdate,format='%Y%m%d%H')
date2 = pd.to_datetime(edate,format='%Y%m%d%H')
delta = timedelta(hours=hint)
dates = pd.date_range(start=date1, end=date2, freq=delta)

rule = rrulewrapper(DAILY, byhour=6, interval=8)
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
        print('Processing Radfile: %s' %(raddfile),flush=1)
        
        if 'chkwvn_list' not in locals():
           ds0=xa.open_dataset(infile0)
           ds0=ds0.swap_dims({"nchans":"wavenumber"})
           chkwvn_list=ds0.wavenumber.sel(wavenumber=spectral_range)[ds0.use_flag.sel(wavenumber=spectral_range)==1]
       
        tmpds0=read_rad_ncdiag(infile0,chkwvn=chkwvn_list)
        tmpds1=read_rad_ncdiag(infile1,chkwvn=chkwvn_list)
        
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

for chkwvn in [962.5,1096]:
#for chkwvn in chkwvn_list:
    qc_chk=usedcnts_all.sel(wavenumber=chkwvn)
    
    fig,ax=plt.subplots()
    set_size(axe_w,axe_h,b=0.2)
    ax.plot_date(dates,qc_chk.exp0cnts,'.-',color='b',linewidth=0.8)
    ax.plot_date(dates,qc_chk.exp1cnts,'.-',color='r',linewidth=0.8)
    
    ax.xaxis.set_major_locator(loc)
    ax.xaxis.set_major_formatter(formatter)
    
    tistr=('%s (%.2f $cm^{-1}$)' %(sensor,chkwvn))
    ax.set_title(tistr,loc='left')
    #ax.set_xlabel('Dates')
    ax.set_ylabel('Number of observations')
    ax.legend(expnlist)
            
    if (fsave):
        fname=('%s/TS_%s_%s_%.2f.%s' %(ts_savedir,area,sensor,chkwvn,ffmt))
        print(fname,flush=1)
        fig.savefig(fname,dpi=quality)
        plt.close()

wvn=usedcnts_all.wavenumber.data
wvl=1e4/wvn
wvnlb='Wavenumber [$cm^{-1}$]' 
wvllb='Wavelength [µm]'
colorlst=['b','r']
lstylelst=[' ',' ']
mrklst=['o','^']
lglst=expnlist
cntsyaxlb='Counts'
aerpyaxlb='Aerosol-affected [%]'
tistr=''
prop_dict={'color'     :['b','r'],
           'line_style':[' ',' '],
           'line_width':[1.5,1.5],
           'marker'    :['o','^'],
           'mark_size' :[5.,5.],
           'fillstyle' :['none','none'],
           'legend'    :expnlist,
           }

used0_mean=usedcnts_all.exp0cnts.mean(dim='dates')
used1_mean=usedcnts_all.exp1cnts.mean(dim='dates')
pltds=xa.Dataset({'exp0':(['channels'],used0_mean.data),
                  'exp1':(['channels'],used1_mean.data),
                  },coords={'channels':chkwvn_list.data})

fig,ax=plt.subplots()
set_size(axe_w,axe_h,ax=ax,b=0.25)
plt_2exps_x2y(pltds,cntsyaxlb,wvn,wvnlb,wvl,wvllb,prop_dict,tistr,0,[],plot_diff=2,stat_type='VALUE')

fname=('%s/Spec_%s_%s_cnts.%s-%s.%s' %(sp_savedir,area,sensor,spectral_range.start,spectral_range.stop,ffmt))
if (fsave):
    print(fname,flush=1)
    fig.savefig(fname,dpi=quality)
    plt.close()

aerp0_mean=usedcnts_all.exp0aerp.mean(dim='dates')*100.
aerp1_mean=usedcnts_all.exp1aerp.mean(dim='dates')*100.
pltds=xa.Dataset({'exp0':(['channels'],aerp0_mean.data),
                  'exp1':(['channels'],aerp1_mean.data),
                  },coords={'channels':chkwvn_list.data})
        
fig,ax=plt.subplots()
set_size(axe_w,axe_h,ax=ax,b=0.25)
plt_2exps_x2y(pltds,aerpyaxlb,wvn,wvnlb,wvl,wvllb,prop_dict,tistr,0,[],stat_type='VALUE')

fname=('%s/Spec_%s_%s_aerp.%s-%s.%s' %(sp_savedir,area,sensor,spectral_range.start,spectral_range.stop,ffmt))
if (fsave):
    print(fname,flush=1)
    fig.savefig(fname,dpi=quality)
    plt.close()
