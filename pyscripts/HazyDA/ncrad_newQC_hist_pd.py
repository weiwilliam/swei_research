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
machine='S4'
os_name=platform.system()
if (machine=='MBP'):
    rootpath='/Users/weiwilliam'
    rootarch='/Volumes/WD2TB/ResearchData'
elif (machine=='Desktop'):
    rootpath='F:\GoogleDrive_NCU\Albany'
    rootarch='F:\ResearchData'
    rootgit='F:\GitHub\swei_research'
elif (machine=='S4'):
    rootarch='/data/users/swei'
    rootpath='/data/users/swei'
    rootgit='/home/swei/research'
elif (machine=='Hera'):
    rootpath='/scratch2/BMC/gsd-fv3-dev/Shih-wei.Wei'
    rootarch='/scratch2/BMC/gsd-fv3-dev/Shih-wei.Wei/ResearchData'
    rootgit='/home/Shih-wei.Wei/research'
sys.path.append(rootgit+'/pyscripts/functions')
import setuparea as setarea
from plot_utils import set_size
from utils import ndate
from datetime import datetime, timedelta
import scipy.stats

tlsize=12 ; lbsize=10
mpl.rc('axes', titlesize=tlsize,labelsize=lbsize)
mpl.rc('xtick',labelsize=lbsize)
mpl.rc('ytick',labelsize=lbsize)
mpl.rc('legend',fontsize='small')
fsave=1 ; ffmt='png' ; ptsize=4
axe_w=3 ; axe_h=3 ; quality=300
minussign=u'\u2212'

# Plotting setup
sdate=2020061000
edate=2020071018
aertype='Dust'
hint=6
exp='aerqc_corR'
sensor='iasi_metop-a'
spectral_range=slice(700,1300)
loop='ges' #ges,anl
#if loop=='anl':
#    tlstr='OMA'
#elif loop=='ges':
#    tlstr='OMF'
lutver='v4'
lutfmt='csv'

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
#archpath=rootarch+'/ResearchData/Prospectus/AeroObsStats/nc_diag'
archpath='/scratch/users/swei/ncdiag'
lutpath=rootpath+'/AlbanyWork/Prospectus/Experiments/AeroObsStats/SD_LUT'
outpath=rootpath+'/AlbanyWork/Prospectus/Experiments/HazyDA/Images'
archdir=archpath+'/'+exp
savedir=outpath+'/aerstat/'+exp+'/hist/'+aertype
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
    raddfile='diag_'+sensor+'_'+loop+'.'+str(date)+'.nc4'
    infile1=archdir+'/'+str(date)+'/'+raddfile

    if (os.path.exists(infile1)):
        print('Processing Radfile: %s' %(raddfile), flush=1)
        ds1=xa.open_dataset(infile1)
        npts=int(ds1.nobs.size/ds1.nchans.size)
        nchs=ds1.nchans.size
        ds1=ds1.swap_dims({"nchans":"wavenumber"}) #replace the dimension of channel by channel indices
        wavelength=1e+04/ds1.wavenumber
        chkwvn_list=ds1.wavenumber.sel(wavenumber=ds1.wavenumber[ds1.use_flag==1])
        dates_count+=1
    else:
        print('%s is not existing'%(raddfile))
        continue
    
    # Observation lat/lon
    rlat1=np.reshape(ds1.Latitude.values,(npts,nchs))[:,0]
    rlon1=np.reshape(ds1.Longitude.values,(npts,nchs))[:,0]
    qcflags=np.reshape(ds1.QC_Flag.values,(npts,nchs))
    #obs1=np.reshape(ds1.Observation.values,(npts,nchs))
    sim1=np.reshape(ds1.Simulated_Tb.values,(npts,nchs))
    clr1=np.reshape(ds1.Clearsky_Tb.values,(npts,nchs))
    varinv1=np.reshape(ds1.Inverse_Observation_Error.values,(npts,nchs))
    omb_nbc1=np.reshape(ds1.Obs_Minus_Forecast_unadjusted.values,(npts,nchs))
    obs1=omb_nbc1+sim1
    aereff_fg=sim1-clr1
    aereff_obs=obs1-clr1
    aereff=0.5*abs(aereff_fg)+0.5*abs(aereff_obs)
    tmpds=xa.Dataset({'rlon1':(['obsloc'],rlon1),
                      'rlat1':(['obsloc'],rlat1),
                      'qcflag':(['obsloc','wavenumber'],qcflags),
                      'Ae':(['obsloc','wavenumber'],aereff),
                      'OMF':(['obsloc','wavenumber'],omb_nbc1)},
                     coords={'obsloc':np.arange(npts),
                             'wavenumber':ds1.wavenumber.values})
    
    tmpds=tmpds.sel(wavenumber=chkwvn_list)
    
    if (date==str(sdate)):
        ds_all=tmpds
    else:
        ds_all=xa.concat((ds_all,tmpds),dim='obsloc')

total_obscounts=ds_all.obsloc.size
ds_all=ds_all.assign_coords(obsloc=np.arange(total_obscounts))

satstats_file=lutpath+'/'+sensor+'_'+str(nchs)+'_stats.'+lutver+'.'+lutfmt
if (lutfmt=='xlsx'):
   lutdf=pd.read_excel(satstats_file)
elif (lutfmt=='csv'):
   lutdf=pd.read_csv(satstats_file)
chk_filter=(lutdf['iuse']==1.0)&(lutdf['Aer_sen']==1.0)
chkwvn_list=lutdf.loc[chk_filter,:]['wavenumber'].tolist()

binsize=0.1
halfbin=0.5*binsize
hist_x_edge=np.arange(-1*halfbin,20.+binsize,binsize)
bin_center=(hist_x_edge+halfbin)[:-1]

#for chkwvn in [962.5]:
for chkwvn in chkwvn_list:
    df_chk=ds_all.sel(wavenumber=chkwvn).to_dataframe()
    tmpdf=lutdf.loc[lutdf['wavenumber']==chkwvn]
    sdmin=tmpdf.SD_min.values[0]
    sdmax=tmpdf.SD_max.values[0]
    ae1=tmpdf.Aeff_1.values[0]
    ae2=tmpdf.Aeff_2.values[0]
    obserr=tmpdf.SD_o.values[0]
    aedepsd=np.zeros_like(bin_center,dtype='float')
    aedepsd[:]=obserr
    aedepsd=np.where(bin_center>=ae2,sdmax,aedepsd)
    for idx in np.where((bin_center<ae2)&(bin_center>ae1))[0]:
        aedepsd[idx]=sdmin+(sdmax-sdmin)/(ae2-ae1)*(bin_center[idx]-ae1)

    filter_ori=(df_chk['qcflag']!=7.0)
    #filter_bst=(abs(df_chk['OMF'])>3)&(abs(df_chk['OMF'])>1.8*df_chk['Ae'])
    filter_fnl=(df_chk['qcflag']==0.0)|(df_chk['qcflag']==13.0)
    
    ori_df=df_chk.loc[filter_ori,:]
    fnl_df=df_chk.loc[filter_fnl,:]
    
    ori_df['Ae_bin']=pd.cut(ori_df['Ae'],bins=hist_x_edge,labels=bin_center)
    fnl_df['Ae_bin']=pd.cut(fnl_df['Ae'],bins=hist_x_edge,labels=bin_center)
    
    ori_bin_df=ori_df.groupby('Ae_bin').agg({'OMF':['count','mean','std']})
    fnl_bin_df=fnl_df.groupby('Ae_bin').agg({'OMF':['count','mean','std']})
    ori_bin_df=ori_bin_df.reset_index()
    fnl_bin_df=fnl_bin_df.reset_index()

    ori_omb_mean=ori_bin_df['OMF','mean']
    fnl_omb_mean=fnl_bin_df['OMF','mean']
    ori_omb_std=ori_bin_df['OMF','std']
    fnl_omb_std=fnl_bin_df['OMF','std']

    x_label='Aerosol effect [K]'
    
    tistr=('%s (%.2f $cm^{-1}$)' %(sensor,chkwvn))
    
    fig=plt.figure()
    ax=plt.subplot()
    set_size(axe_w,axe_h,l=0.15,r=0.85)
    ax.set_title(tistr,loc='left')
    ax.set_xlabel(x_label)
    ax.bar(bin_center,ori_bin_df['OMF','count'],binsize,color='grey',alpha=0.3)
    ax.bar(bin_center,fnl_bin_df['OMF','count'],binsize,color='grey',alpha=0.7)
    ax.set_yscale("log")
    ax.set_ylabel("Counts")
    ax2=ax.twinx()
    ax2.plot(bin_center,ori_omb_mean,'tab:blue',linewidth=0.7,label='Mean')
    ax2.plot(bin_center,ori_omb_std,'tab:gray',linewidth=0.7,label='SD')
    ax2.plot(bin_center,fnl_omb_mean,'tab:blue',label='Mean_QC')
    ax2.plot(bin_center,fnl_omb_std,'tab:red',label='SD_QC')
    ax2.plot(bin_center,aedepsd,'k',linewidth=1.,label='${A}_{e}$_dep_SD')
    ax2.set_ylim(-40,10)
    ax2.set_ylabel('SD and Mean of O%sF [K]'%(minussign))
    ax2.hlines(obserr,0,1,transform=ax2.get_yaxis_transform(),colors='k',linestyle='dashed',linewidth=0.7)
    ax2.hlines(0.,0,1,transform=ax2.get_yaxis_transform(),colors='grey',linewidth=0.4)
    ax2.legend(loc=7)
#    
    if (fsave):
        fname=('%s/PDF_MeanSD_%s_%s_%.2f.%s'
                %(savedir,area,sensor,chkwvn,ffmt))
        print(fname,flush=1)
        fig.savefig(fname,dpi=quality)
        plt.close()
