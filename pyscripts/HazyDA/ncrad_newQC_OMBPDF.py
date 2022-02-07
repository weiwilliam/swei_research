# -*- coding: utf-8 -*-
"""
Created on Tue Jun 29 16:06:23 2021

@author: ck102

Aerosol detection based on CADS 3.1 from NWP SAF

"""
import sys, os, platform
import numpy as np
import xarray as xa
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mpcrs
import cartopy.crs as ccrs
os_name=platform.system()
if (os_name=='Darwin'):
    rootpath='/Users/weiwilliam'
    rootarch='/Volumes/WD2TB/ResearchData'
elif (os_name=='Windows'):
    rootpath='F:\GoogleDrive_NCU\Albany'
    rootarch='F:\ResearchData'
    rootgit='F:\GitHub\swei_research'
sys.path.append(rootgit+'/pyscripts/functions')
import setuparea as setarea
from plot_utils import setupax_2dmap, plt_x2y, set_size
from utils import ndate,setup_cmap
from datetime import datetime, timedelta
import scipy.stats

tlsize=12 ; lbsize=10
mpl.rc('axes', titlesize=tlsize,labelsize=lbsize)
mpl.rc('xtick',labelsize=lbsize)
mpl.rc('ytick',labelsize=lbsize)
mpl.rc('legend',fontsize='small')
fsave=1 ; ffmt='png' ; ptsize=4
axe_w=3 ; axe_h=3 ; quality=300

# Projection setting
proj=ccrs.PlateCarree(globe=None)

# Plotting setup
sdate=2020082200
edate=2020092118
aertype='Smoke'
hint=6
exp0='GDAS'
exp1='aerqc_corR'
lglst=['Clear','Hazy']
sensor='iasi_metop-a'
spectral_range=slice(600,1300)
loop='ges' #ges,anl
usebc=0
plthist=1 # plot 2d histogram

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
   
if (usebc):
    bcflg='bc'
else:
    bcflg='nobc'

# Data path setup
archpath=rootarch+'/Prospectus/AeroObsStats/nc_diag'
lutpath='F:\GoogleDrive_NCU\Albany\AlbanyWork\Prospectus\Experiments\AeroObsStats\SD_LUT\All'
outpath=rootpath+'/AlbanyWork/Prospectus/Experiments/AeroObsStats/images/'
archdir0=archpath+'/'+exp0
archdir1=archpath+'/'+exp1
savedir=outpath+'/'+exp1+'/pdf/'+aertype
if ( not os.path.exists(savedir) ):
    os.makedirs(savedir)

syy=int(str(sdate)[:4]); smm=int(str(sdate)[4:6])
sdd=int(str(sdate)[6:8]); shh=int(str(sdate)[8:10])
eyy=int(str(edate)[:4]); emm=int(str(edate)[4:6])
edd=int(str(edate)[6:8]); ehh=int(str(edate)[8:10])

date1 = datetime(syy,smm,sdd,shh)
date2 = datetime(eyy,emm,edd,ehh)
delta = timedelta(hours=hint)
dates = pd.date_range(start=date1, end=date2, freq=delta)

tnum=0
dlist=[]
cdate=sdate
while (cdate<=edate):
    dlist.append(str(cdate))
    tnum=tnum+1
    cdate=ndate(hint,cdate)

dates_count=0
for date in dlist:
    raddfile0='diag_'+sensor+'_'+loop+'.'+str(date)+'.nc'
    raddfile1='diag_'+sensor+'_'+loop+'.'+str(date)+'.nc4'
    infile0=archdir0+'/'+raddfile0
    infile1=archdir1+'/'+str(date)+'/'+raddfile1

    if (os.path.exists(infile0) and 
        os.path.exists(infile1)):
        print('Processing Radfile: %s' %(raddfile1))
        ds0=xa.open_dataset(infile0)
        ds1=xa.open_dataset(infile1)
        npts0=ds0.obsloc.size
        npts1=int(ds1.nobs.size/ds1.nchans.size)
        nchs0=ds0.channel.size
        nchs1=ds1.nchans.size
        ds0=ds0.assign_coords(nuchan=('wavenumber',ds0.wavenumber))
        ds0=ds0.swap_dims({"channel":"wavenumber"})
        ds1=ds1.assign_coords(nuchan=('wavenumber',ds1.wavenumber))
        ds1=ds1.swap_dims({"nchans":"wavenumber"}) #replace the dimension of channel by channel indices
        wavelength=1e+04/ds1.wavenumber
        chkwvn_list=ds1.wavenumber.sel(wavenumber=spectral_range)[ds1.use_flag.sel(wavenumber=spectral_range)==1]
        #usedchidx=np.where(ds1.iuse_rad==1)[0]
        #unusedchidx=np.where(ds1.iuse_rad==-1)[0]
        #wvldiff=abs(np.subtract(wavelength[usedchidx],chkwvl))
        #chkwvlidx=usedchidx[np.where(wvldiff==wvldiff.min())[0][0]]
        #print('Check wavelength: %.2f' %(wavelength[chkwvlidx]))
        dates_count+=1
    else:
        print('%s is not existing'%(raddfile1))
        continue

    # Observation lat/lon from exp 0 (baseline)
    if (exp0=='GDAS'):
        rlat0=ds0.locinfo[0,:]
        rlon0=ds0.locinfo[1,:]
        qcflags0=ds0.qcflag
        obs0=ds0.tb_obs
        sim0=obs0-ds0.tbcnob
        clr0=xa.zeros_like(sim0)
        omb_bc0=ds0.tbc
        omb_nbc0=ds0.tbcnob
        varinv0=ds0.errinv
        arrshape=['wavenumber','obsloc']
    else:
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
    tmpds0=tmpds0.sel(wavenumber=spectral_range)
   
      
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
    tmpds1=tmpds1.sel(wavenumber=spectral_range)
    
    if (date==str(sdate)):
        ds_all0=tmpds0
        ds_all1=tmpds1
    else:
        ds_all0=xa.concat((ds_all0,tmpds0),dim='obsloc')
        ds_all1=xa.concat((ds_all1,tmpds1),dim='obsloc')

total_obscounts=ds_all1.obsloc.size
ds_all1=ds_all1.assign_coords(obsloc=np.arange(total_obscounts))

# satinfo_excel=lutpath+'/'+sensor+'_'+str(nchs)+'_stats.xlsx'
# lutdf=pd.read_excel(satinfo_excel)
satinfo_csv=lutpath+'/'+sensor+'_'+str(nchs1)+'_stats_new.v3.csv'
lutdf=pd.read_csv(satinfo_csv)
filter = ((lutdf.Aer_sen==1.)&(lutdf.iuse==1.)&
          (lutdf.wavenumber>=spectral_range.start)&
          (lutdf.wavenumber<=spectral_range.stop))
tmpdf=lutdf.loc[filter,:]
chkwvn_list=tmpdf.wavenumber.values

binsize=0.1
halfbin=0.5*binsize
hist_x_edge=np.arange(-10-halfbin,10.+binsize,binsize)
bin_center=(hist_x_edge+halfbin)[:-1]

# for chkwvn in [962.5]:
for chkwvn in chkwvn_list:
    ds_chk0=ds_all0.sel(wavenumber=chkwvn)
    tb_obs0=ds_chk0.tb_obs
    # varinv1=ds_chk1.varinv
    
    ds_chk1=ds_all1.sel(wavenumber=chkwvn)
    tb_obs1=ds_chk1.tb_obs
    # varinv1=ds_chk1.varinv
    
    if (usebc):
        omb0=ds_chk0.omb_bc
        omb1=ds_chk1.omb_bc
    else:
        omb0=ds_chk0.omb_nbc
        omb1=ds_chk1.omb_nbc
    # aereff_fg=tb_sim-tb_clr
    # aereff_obs=tb_obs-tb_clr
    # aereff=0.5*abs(aereff_fg)+0.5*abs(aereff_obs)

    # tmpdf=lutdf.loc[lutdf['wavenumber']==chkwvn]
    # sdmin=tmpdf.SD_min.values[0]
    # sdmax=tmpdf.SD_max.values[0]
    # ae1=tmpdf.Aeff_1.values[0]
    # ae2=tmpdf.Aeff_2.values[0]
    # obserr=tmpdf.SD_o.values[0]
    
    good_msk0=(ds_chk0.qcflag==0.)
    gross_msk0=(ds_chk0.qcflag==3.)
    cld_msk0=(ds_chk0.qcflag==7.)
    tzr_msk0=(ds_chk0.qcflag==10.)
    sfcir_msk0=(ds_chk0.qcflag==53.)
    passed_msk0=(good_msk0)
    # ori_msk=(good_msk)|(aer_msk)|(bust_msk)|(tzr_msk)

    good_msk1=(ds_chk1.qcflag==0.)
    gross_msk1=(ds_chk1.qcflag==3.)
    cld_msk1=(ds_chk1.qcflag==7.)
    tzr_msk1=(ds_chk1.qcflag==10.)
    aer_msk1=(ds_chk1.qcflag==13.)
    sfcir_msk1=(ds_chk1.qcflag==53.)
    bust_msk1=(ds_chk1.qcflag==55.)
    passed_msk1=(good_msk1)|(aer_msk1)
    # ori_msk=(good_msk)|(aer_msk)|(bust_msk)|(tzr_msk)
    
    pltmsk0=passed_msk0
    pltmsk1=passed_msk1

    if (plthist):            
        pltda_x0=omb0[pltmsk0==1]
        pltda_x1=omb1[pltmsk1==1]
        
        x_label='OMB [K]'
        hdata0, tmpbins=np.histogram(pltda_x0, bins=hist_x_edge, density=1)
        hdata1, tmpbins=np.histogram(pltda_x1, bins=hist_x_edge, density=1)
    
        omb_mean=np.zeros_like(bin_center,dtype='float')
        omb_sd=np.zeros_like(bin_center,dtype='float')
        # counts=np.zeros_like(bin_center,dtype='int')
        
        tistr=('%s (%.2f $cm^{-1}$)' %(sensor,chkwvn))
        
        fig=plt.figure()
        ax=plt.subplot()
        set_size(axe_w,axe_h,l=0.15,r=0.85)
        ax.set_title(tistr,loc='left')
        ax.set_xlabel(x_label)
        ax.plot(bin_center,hdata0*binsize,'tab:blue')
        ax.plot(bin_center,hdata1*binsize,'tab:red')
        # ax.plot(normbin,normdis,'k--')
        ax.set_yscale("log")
        ax.set_ylabel("PDF")
        ax.set_ylim(1e-4,1.4e0)
        ax.vlines(0.,0,1,transform=ax.get_xaxis_transform(),colors='grey',linestyle='dashed',linewidth=0.7)
        ax.legend(lglst,loc=2)
        
        if (fsave):
            fname=('PDF_OMB_%s_%s_%s.%.2f.%s'
                    %(area,sensor,bcflg,chkwvn,ffmt))
            fig.savefig(savedir+'/'+fname,dpi=quality)
            plt.close()
