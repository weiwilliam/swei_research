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
    rootpath='/data/users/swei/'
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

# Plotting setup
sdate=2020061000
edate=2020071018
hint=6
explist=['hazyda_ctrl','hazyda_aero']
expnlist=['CTL','AER']
sensor='iasi_metop-a'
spectral_range=slice(600,1300)
loop='ges' #ges,anl
#if loop=='anl':
#    tlstr='OMA'
#elif loop=='ges':
#    tlstr='OMF'
plthist=1 # plot 2d histogram
usebc=0

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
outpath=rootpath+'/AlbanyWork/Prospectus/Experiments/HazyDA/Images'
archdir0=rootarch+'/'+explist[0]
archdir1=rootarch+'/'+explist[1]
savedir=outpath+'/DiagFiles/rad/pdf_norm/'+area
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
    raddfile0='diag_'+sensor+'_'+loop+'.'+str(date)+'.nc4'
    raddfile1='diag_'+sensor+'_'+loop+'.'+str(date)+'.nc4'
    infile0=archdir0+'/'+str(date)+'/'+raddfile0
    infile1=archdir1+'/'+str(date)+'/'+raddfile1

    if (os.path.exists(infile0) and 
        os.path.exists(infile1) ):
        print('Processing Radfile: %s' %(raddfile0),flush=1)
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
        chkwvn_list=ds1.wavenumber[ds1.use_flag==1]
        dates_count+=1
    else:
        print('%s is not existing'%(raddfile0),flush=1)
        continue
    
    # Observation lat/lon from exp 0 (baseline)
    rlat0=np.reshape(ds0.Latitude.values,(npts0,nchs0))
    rlon0=np.reshape(ds0.Longitude.values,(npts0,nchs0))
    qcflags0=np.reshape(ds0.QC_Flag.values,(npts0,nchs0))
    #obs0=np.reshape(ds0.Observation.values,(npts0,nchs0))
    sim0=np.reshape(ds0.Simulated_Tb.values,(npts0,nchs0))
    clr0=np.reshape(ds0.Clearsky_Tb.values,(npts0,nchs0))
    varinv0=np.reshape(ds0.Inverse_Observation_Error.values,(npts0,nchs0))
    sim_bc0=np.reshape(ds0.Obs_Minus_Forecast_adjusted.values,(npts0,nchs0))
    sim_nbc0=np.reshape(ds0.Obs_Minus_Forecast_unadjusted.values,(npts0,nchs0))
    obs0=sim_nbc0+sim0
    tmpds0=xa.Dataset({'rlon':(['obsloc'],rlon0[:,0]),
                      'rlat':(['obsloc'],rlat0[:,0]),
                      'qcflag':(['obsloc','wavenumber'],qcflags0),
                      'tb_obs':(['obsloc','wavenumber'],obs0),
                      'tb_sim':(['obsloc','wavenumber'],sim0),
                      'tb_clr':(['obsloc','wavenumber'],clr0),
                      'varinv':(['obsloc','wavenumber'],varinv0),
                      'omb_bc':(['obsloc','wavenumber'],sim_bc0),
                      'omb_nbc':(['obsloc','wavenumber'],sim_nbc0)},
                      coords={'obsloc':np.arange(npts0),
                             'wavenumber':ds0.wavenumber.values})
    tmpds0=tmpds0.sel(wavenumber=chkwvn_list)
    
    # Observation lat/lon from exp 1 (test)
    rlat1=np.reshape(ds1.Latitude.values,(npts1,nchs1))
    rlon1=np.reshape(ds1.Longitude.values,(npts1,nchs1))
    qcflags1=np.reshape(ds1.QC_Flag.values,(npts1,nchs1))
    #obs1=np.reshape(ds1.Observation.values,(npts1,nchs1))
    sim1=np.reshape(ds1.Simulated_Tb.values,(npts1,nchs1))
    clr1=np.reshape(ds1.Clearsky_Tb.values,(npts1,nchs1))
    varinv1=np.reshape(ds1.Inverse_Observation_Error.values,(npts1,nchs1))
    sim_bc1=np.reshape(ds1.Obs_Minus_Forecast_adjusted.values,(npts1,nchs1))
    sim_nbc1=np.reshape(ds1.Obs_Minus_Forecast_unadjusted.values,(npts1,nchs1))
    obs1=sim_nbc1+sim1
    tmpds1=xa.Dataset({'rlon':(['obsloc'],rlon1[:,0]),
                      'rlat':(['obsloc'],rlat1[:,0]),
                      'qcflag':(['obsloc','wavenumber'],qcflags1),
                       'tb_obs':(['obsloc','wavenumber'],obs1),
                       'tb_sim':(['obsloc','wavenumber'],sim1),
                       'tb_clr':(['obsloc','wavenumber'],clr1),
                      'varinv':(['obsloc','wavenumber'],varinv1),
                      'omb_bc':(['obsloc','wavenumber'],sim_bc1),
                      'omb_nbc':(['obsloc','wavenumber'],sim_nbc1)},
                      coords={'obsloc':np.arange(npts1),
                             'wavenumber':ds1.wavenumber.values})
    tmpds1=tmpds1.sel(wavenumber=chkwvn_list)
    
    if (date==str(sdate)):
        ds_all0=tmpds0
        ds_all1=tmpds1
    else:
        ds_all0=xa.concat((ds_all0,tmpds0),dim='obsloc')
        ds_all1=xa.concat((ds_all1,tmpds1),dim='obsloc')

total_obscounts=ds_all0.obsloc.size
ds_all0=ds_all0.assign_coords(obsloc=np.arange(total_obscounts))
ds_all1=ds_all1.assign_coords(obsloc=np.arange(total_obscounts))

binsize=0.1
halfbin=0.5*binsize
hist_x_edge=np.arange(-3-halfbin,3.+binsize,binsize)
bin_center=(hist_x_edge+halfbin)[:-1]

# for chkwvn in [652.0]:
for chkwvn in [962.5]:
# for chkwvn in chkwvn_list:
    ds_chk0=ds_all0.sel(wavenumber=chkwvn)
    tb_sim0=ds_chk0.tb_sim
    tb_clr0=ds_chk0.tb_clr
    tb_obs0=ds_chk0.tb_obs
    omb_bc0=ds_chk0.omb_bc
    omb_nbc0=ds_chk0.omb_nbc
    varinv0=ds_chk0.varinv
    qcflg0=ds_chk0.qcflag

    ds_chk1=ds_all1.sel(wavenumber=chkwvn)
    tb_sim1=ds_chk1.tb_sim
    tb_clr1=ds_chk1.tb_clr
    tb_obs1=ds_chk1.tb_obs
    omb_bc1=ds_chk1.omb_bc
    omb_nbc1=ds_chk1.omb_nbc
    varinv1=ds_chk1.varinv
    qcflg1=ds_chk1.qcflag
    
    omb0=omb_nbc0
    omb1=omb_nbc1
    
    good_msk0=(ds_chk0.qcflag==0.)
    gross_msk0=(ds_chk0.qcflag==3.)
    cld_msk0=(ds_chk0.qcflag==7.)
    tzr_msk0=(ds_chk0.qcflag==10.)
    aer_msk0=(ds_chk0.qcflag==13.)
    sfcir_msk0=(ds_chk0.qcflag==53.)
    bust_msk0=(ds_chk0.qcflag==55.)
    passed_msk0=(good_msk0)|(aer_msk0)
    ori_msk0=(good_msk0)|(aer_msk0)|(bust_msk0)|(tzr_msk0)|(gross_msk0)
    pltmsk0=passed_msk0
    
    good_msk1=(ds_chk1.qcflag==0.)
    gross_msk1=(ds_chk1.qcflag==3.)
    cld_msk1=(ds_chk1.qcflag==7.)
    tzr_msk1=(ds_chk1.qcflag==10.)
    aer_msk1=(ds_chk1.qcflag==13.)
    sfcir_msk1=(ds_chk1.qcflag==53.)
    bust_msk1=(ds_chk1.qcflag==55.)
    passed_msk1=(good_msk1)|(aer_msk1)
    ori_msk1=(good_msk1)|(aer_msk1)|(bust_msk1)|(tzr_msk1)|(gross_msk1)
    pltmsk1=passed_msk1

    if (plthist):            
        pltda_x1=omb0[pltmsk0==1]*varinv0[pltmsk0==1]
        pltda_x2=omb1[pltmsk1==1]*varinv1[pltmsk1==1]
        # pltda_x4=omb[final_qc_msk==1]/obs_sd2[final_qc_msk==1]
        x_label='Normalized OMB'
        hdata1, tmpbins=np.histogram(pltda_x1, bins=hist_x_edge, density=0)
        hdata2, tmpbins=np.histogram(pltda_x2, bins=hist_x_edge, density=0)
    
        omb_mean=np.zeros_like(bin_center,dtype='float')
        omb_sd=np.zeros_like(bin_center,dtype='float')
        # counts=np.zeros_like(bin_center,dtype='int')
        
        tistr=('%s (%.2f $cm^{-1}$)' %(sensor,chkwvn))
        
        fig=plt.figure()
        ax=plt.subplot()
        set_size(axe_w,axe_h,l=0.15,r=0.85)
        ax.set_title(tistr,loc='left')
        ax.set_xlabel(x_label)
        ax.plot(bin_center,hdata1*binsize,'b')
        ax.plot(bin_center,hdata2*binsize,'r')
        #ax.set_yscale("log")
        #ax.set_ylim(1e-4,1.4e0)
        #ax.set_ylabel("PDF")
        ax.set_ylabel('Counts')
        ax.vlines(0.,0,1,transform=ax.get_xaxis_transform(),colors='grey',linestyle='dashed',linewidth=0.7)
        ax.legend(expnlist,loc=2)
        
        if (fsave):
            fname=('PDF_NormOMB_%s_%s_%.2f.%s_%s.%s'
                    %(area,sensor,chkwvn,sdate,edate,ffmt))
            print(fname,flush=1)
            fig.savefig(savedir+'/'+fname,dpi=quality)
            plt.close()
