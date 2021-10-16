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
sys.path.append(rootpath+'/AlbanyWork/Utility/Python3/functions')
import setuparea as setarea
from plot_utils import setupax_2dmap, plt_x2y, set_size
from utils import ndate,setup_cmap
from datetime import datetime, timedelta
import scipy.stats

tlsize=12 ; lbsize=10
mpl.rc('axes', titlesize=tlsize,labelsize=lbsize)
mpl.rc('xtick',labelsize=lbsize)
mpl.rc('ytick',labelsize=lbsize)
mpl.rc('legend',fontsize='large')
fsave=1 ; ffmt='png' ; ptsize=4
axe_w=3 ; axe_h=3 ; quality=300

# Projection setting
proj=ccrs.PlateCarree(globe=None)

# Plotting setup
sdate=2020061000
edate=2020071018
aertype='Dust'
hint=6
exp='AerObserver'
sensor='iasi_metop-a'
spectral_range=slice(700,1300)
loop='ges' #ges,anl
#if loop=='anl':
#    tlstr='OMA'
#elif loop=='ges':
#    tlstr='OMF'
plthist=0 # plot 2d histogram
plthist_mean_sd=1 # plot 2d histogram

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
archpath=rootarch+'/Prospectus/AeroObsStats/nc_diag'
lutpath='F:\GoogleDrive_NCU\Albany\AlbanyWork\Prospectus\Experiments\AeroObsStats\images\AerObserver_newQC\SD_LUT\All'
outpath=rootpath+'/AlbanyWork/Prospectus/Experiments/AeroObsStats/images/'
archdir=archpath+'/'+exp
savedir=outpath+'/'+exp+'_newQC/pdf/'+aertype
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
        print('Processing Radfile: %s' %(raddfile))
        ds1=xa.open_dataset(infile1)
        npts=int(ds1.nobs.size/ds1.nchans.size)
        nchs=ds1.nchans.size
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
        print('%s is not existing'%(raddfile))
        continue
    
    # Observation lat/lon
    rlat1=np.reshape(ds1.Latitude.values,(npts,nchs))
    rlon1=np.reshape(ds1.Longitude.values,(npts,nchs))
    qcflags=np.reshape(ds1.QC_Flag.values,(npts,nchs))
    obs1=np.reshape(ds1.Observation.values,(npts,nchs))
    sim1=np.reshape(ds1.Simulated_Tb.values,(npts,nchs))
    clr1=np.reshape(ds1.Clearsky_Tb.values,(npts,nchs))
    varinv1=np.reshape(ds1.Inverse_Observation_Error.values,(npts,nchs))
    tmpds=xa.Dataset({'rlon1':(['obsloc'],rlon1[:,0]),
                      'rlat1':(['obsloc'],rlat1[:,0]),
                      'qcflag':(['obsloc','wavenumber'],qcflags),
                      'tb_obs':(['obsloc','wavenumber'],obs1),
                      'tb_sim':(['obsloc','wavenumber'],sim1),
                      'tb_clr':(['obsloc','wavenumber'],clr1),
                      'varinv':(['obsloc','wavenumber'],varinv1)},
                     coords={'obsloc':np.arange(npts),
                             'wavenumber':ds1.wavenumber.values})
    
    tmpds=tmpds.sel(wavenumber=spectral_range)
    
    if (date==str(sdate)):
        ds_all=tmpds
    else:
        ds_all=xa.concat((ds_all,tmpds),dim='obsloc')

total_obscounts=ds_all.obsloc.size
ds_all=ds_all.assign_coords(obsloc=np.arange(total_obscounts))

satinfo_excel=lutpath+'/'+sensor+'_'+str(nchs)+'_stats.xlsx'
lutdf=pd.read_excel(satinfo_excel)

binsize=0.1
halfbin=0.5*binsize
hist_x_edge=np.arange(-3,3.1,binsize)
bin_center=(hist_x_edge+binsize*0.5)[:-1]

normbin=np.arange(-3,3.1,0.1)
normdis_mean=0.
normdis_std=0.5
normdis=scipy.stats.norm.pdf(normbin,normdis_mean,normdis_std)

# for chkwvn in [962.5]:
for chkwvn in chkwvn_list.values:
    ds_chk=ds_all.sel(wavenumber=chkwvn)
    tb_sim=ds_chk.tb_sim
    tb_clr=ds_chk.tb_clr
    tb_obs=ds_chk.tb_obs
    varinv=ds_chk.varinv

    omb=tb_obs-tb_sim
    aereff_fg=tb_sim-tb_clr
    aereff_obs=tb_obs-tb_clr
    aereff=0.5*abs(aereff_fg)+0.5*abs(aereff_obs)

    tmpdf=lutdf.loc[lutdf['wavenumber']==chkwvn]
    sdmin=tmpdf.SD_min.values[0]
    sdmax=tmpdf.SD_max.values[0]
    ae1=tmpdf.Aeff_1.values[0]
    ae2=tmpdf.Aeff_2.values[0]
    obserr=tmpdf.SD_o.values[0]
    
    qc0_msk=(ds_chk.qcflag==0.)
    # qc7_msk=(ds_chk.qcflag==7.)
    qc13_msk=(ds_chk.qcflag==13.)

    ori_msk=((qc0_msk)|(qc13_msk))
    ori_total=np.count_nonzero(ori_msk)
    # lowbt_qc=(used_msk)&(abs(omb)<30.)
    final_qc_msk=(ori_msk)&(~((abs(omb)>3)&(abs(omb)>1.8*aereff)))&(abs(omb)<30.)

    obs_sd1=xa.zeros_like(omb,dtype='float')
    obs_sd2=xa.zeros_like(omb,dtype='float')
    obs_sd1[:]=omb[ori_msk==1].std()
    le_ae1_msk=(aereff<=ae1)
    betweenmsk=((aereff>ae1)&(aereff<ae2))
    ge_ae2_msk=(aereff>=ae2)
    obs_sd2[le_ae1_msk]=sdmin
    obs_sd2[ge_ae2_msk]=sdmax
    for obsidx in ds_chk.obsloc[betweenmsk].values:
        obs_sd2[obsidx]=sdmin+(sdmax-sdmin)/(ae2-ae1)*(aereff[obsidx].values-ae1)

    if (plthist):
        savedir=outpath+'/'+exp+'_newQC/pdf/'+aertype
        if ( not os.path.exists(savedir) ):
            os.makedirs(savedir)
            
        # pltda_x1=omb[ori_msk==1]*varinv[ori_msk==1]
        # pltda_x2=omb[final_qc_msk==1]*varinv[final_qc_msk==1]
        # pltda_x1=omb[ori_msk==1]/obs_sd1[ori_msk==1]
        pltda_x1=omb[ori_msk==1]/omb[ori_msk==1].std()
        pltda_x2=omb[final_qc_msk==1]/obs_sd2[final_qc_msk==1]
        x_label='Normalized OMB'
        hdata1, tmpbins=np.histogram(pltda_x1, bins=hist_x_edge, density=1)
        hdata2, tmpbins=np.histogram(pltda_x2, bins=hist_x_edge, density=1)
    
        omb_mean=np.zeros_like(bin_center,dtype='float')
        omb_sd=np.zeros_like(bin_center,dtype='float')
        counts=np.zeros_like(bin_center,dtype='int')
        
        tistr=('%s (%.2f $cm^{-1}$)' %(sensor,chkwvn))
        
        fig=plt.figure()
        ax=plt.subplot()
        set_size(axe_w,axe_h,l=0.15,r=0.85)
        ax.set_title(tistr,loc='left')
        ax.set_xlabel(x_label)
        ax.plot(bin_center,hdata1*binsize,'tab:blue')
        ax.plot(bin_center,hdata2*binsize,'tab:red')
        ax.plot(normbin,normdis,'k--')
        ax.set_yscale("log")
        ax.set_ylabel("PDF")
        ax.set_ylim(1e-4,1.4e0)
        ax.vlines(0.,0,1,transform=ax.get_xaxis_transform(),colors='grey',linestyle='dashed',linewidth=0.7)
        
        if (fsave):
            fname=('PDF_NormOMB_%s_%s_%.2f.%s'
                    %(area,sensor,chkwvn,ffmt))
            fig.savefig(savedir+'/'+fname,dpi=quality)
            plt.close()

    if (plthist_mean_sd):
        savedir=outpath+'/'+exp+'_newQC/hist/'+aertype
        if ( not os.path.exists(savedir) ):
            os.makedirs(savedir)
        
        pltda_x1=aereff[ori_msk==1]
        pltda_x2=aereff[final_qc_msk==1]
        x_label='Aerosol effect [K]'
        hist_x_edge=np.arange(-1*halfbin,50.+binsize,binsize)
        bin_center=(hist_x_edge+halfbin)[:-1]
        
        omb_mean1=np.zeros_like(bin_center,dtype='float')
        omb_sd1=np.zeros_like(bin_center,dtype='float')
        counts1=np.zeros_like(bin_center,dtype='int')
        omb_mean2=np.zeros_like(bin_center,dtype='float')
        omb_sd2=np.zeros_like(bin_center,dtype='float')
        counts2=np.zeros_like(bin_center,dtype='int')
        
        for i in np.arange(omb_mean1.size):
            lb_aereff=hist_x_edge[i]
            ub_aereff=hist_x_edge[i+1]
            tmpmsk1=(ori_msk)&((aereff>=lb_aereff)&(aereff<ub_aereff))
            omb_mean1[i]=omb[tmpmsk1==1].mean()
            omb_sd1[i]=omb[tmpmsk1==1].std()
            counts1[i]=np.count_nonzero(tmpmsk1)
            tmpmsk2=(final_qc_msk)&((aereff>=lb_aereff)&(aereff<ub_aereff))
            omb_mean2[i]=omb[tmpmsk2==1].mean()
            omb_sd2[i]=omb[tmpmsk2==1].std()
            counts2[i]=np.count_nonzero(tmpmsk2)

        tistr=('%s (%.2f $cm^{-1}$)' %(sensor,chkwvn))
        
        fig=plt.figure()
        ax=plt.subplot()
        set_size(axe_w,axe_h,l=0.15,r=0.85)
        ax.set_title(tistr,loc='left')
        ax.set_xlabel(x_label)
        hdata, bins, patches=ax.hist(pltda_x1,hist_x_edge,color='grey',
                                      density=0,alpha=0.3)
        hdata, bins, patches=ax.hist(pltda_x2,hist_x_edge,color='grey',
                                      density=0,alpha=0.7)
        ax.set_yscale("log")
        ax.set_ylabel("Counts")
        ax2=ax.twinx()
        ax2.plot(bin_center,omb_mean1,'tab:blue',linewidth=0.7)
        ax2.plot(bin_center,omb_sd1,'tab:red',linewidth=0.7)
        ax2.plot(bin_center,omb_mean2,'tab:blue')
        ax2.plot(bin_center,omb_sd2,'tab:red')
        ax2.set_ylim(-40,10)
        ax2.set_ylabel('SD and Mean of O-B [K]')
        ax2.hlines(obserr,0,1,transform=ax2.get_yaxis_transform(),colors='k',linestyle='dashed',linewidth=0.7)
        ax2.hlines(0.,0,1,transform=ax2.get_yaxis_transform(),colors='grey',linewidth=0.4)
        
        if (fsave):
            fname=('PDF_MeanSD_%s_%s_%.2f.%s'
                    %(area,sensor,chkwvn,ffmt))
            fig.savefig(savedir+'/'+fname,dpi=quality)
            plt.close()