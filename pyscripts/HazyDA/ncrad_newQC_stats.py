# -*- coding: utf-8 -*-
"""
Created on Tue Jun 29 16:06:23 2021

@author: ck102

Aerosol detection based on CADS 3.1 from NWP SAF

"""
import sys, os, platform
machine='S4'
import numpy as np
import xarray as xa
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mpcrs
import cartopy.crs as ccrs
os_name=platform.system()
if (machine=='MBP'):
    rootpath='/Users/weiwilliam'
    rootarch='/Volumes/WD2TB/ResearchData'
elif (machine=='Desktop'):
    rootpath='F:\GoogleDrive_NCU\Albany'
    rootarch='F:\ResearchData'
    rootgit='F:\GitHub\swei_research'
elif (machine=='S4'):
    rootarch='/data/users/swei/Experiments'
    rootpath='/data/users/swei'
    rootgit='/home/swei/research'
sys.path.append(rootgit+'/pyscripts/functions')
import setuparea as setarea
from plot_utils import setupax_2dmap, plt_x2y, set_size
from utils import ndate,setup_cmap
from datetime import datetime, timedelta
import cProfile, pstats, io
from pstats import SortKey

#pr = cProfile.Profile()
#pr.enable()

tlsize=12 ; lbsize=10
mpl.rc('axes', titlesize=tlsize,labelsize=lbsize)
mpl.rc('xtick',labelsize=lbsize)
mpl.rc('ytick',labelsize=lbsize)
mpl.rc('legend',fontsize='large')
fsave=1 ; ffmt='png' ; ptsize=4
axe_w=3 ; axe_h=3 ; quality=300

# Plotting setup
sdate=2020061000
edate=2020092118
aertype='All'
hint=6
exp='aerqc_corR'
sensor='iasi_metop-a'
spectral_range=slice(700,1300)
loop='ges' #ges,anl
binsize=0.1
version='v4'
n_seg=8
seg=7

# Data path setup
archpath=rootarch+'/AeroObsStats/OUTPUT'
fixpath=rootgit+'/GSI_exps/fix'
outpath=rootarch+'/AeroObsStats'
archdir=archpath+'/'+exp
print(archpath)
print(fixpath)
print(outpath)
print(archdir)
savedir=outpath+'/SD_LUT'
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

satinfof=fixpath+'/'+sensor+'_satinfo.txt'
fsat=open(satinfof)
errlst=[]
nuch=[]
iuse=[]
for line in fsat.readlines():
    nuch.append(line.split()[1])
    iuse.append(line.split()[2])
    errlst.append(line.split()[3])
fsat.close()
err_array=np.array(errlst,dtype='float')
iuse_array=np.array(iuse,dtype='int')
nuch_array=np.array(nuch,dtype='int')

date=dlist[0]
raddfile='diag_'+sensor+'_'+loop+'.'+str(date)+'.nc4'
infile1=archdir+'/'+str(date)+'/'+raddfile
if (os.path.exists(infile1)):
    print('Processing the first Radfile: %s' %(raddfile),flush=1)
    ds1=xa.open_dataset(infile1)
    npts=int(ds1.nobs.size/ds1.nchans.size)
    nchs=ds1.nchans.size
    ds1=ds1.swap_dims({"nchans":"wavenumber"}) #replace the dimension of channel by channel indices
    wavelength=1e+04/ds1.wavenumber
    chkwvn_list=ds1.wavenumber.values
    ds_sensorinfo=xa.Dataset({'obserr':(['wavenumber'],err_array),
                              'nuchan':(['wavenumber'],nuch_array),
                              'iuse':(['wavenumber'],iuse_array)},
                             coords={'wavenumber':ds1.wavenumber.values})
else:
    print('%s is not existing'%(raddfile),flush=1)

l_seg=int(nchs/n_seg)
seg_ptarr=np.arange(0,nchs,l_seg)
seg_s=seg_ptarr[seg]
if (seg+1==n_seg):
   seg_e=nchs
else:
   seg_e=seg_ptarr[seg+1]

icount=0
for chkwvn in chkwvn_list[seg_s:seg_e]:
#for chkwvn in [962.5]:
    print('Processing for channel with wavenumber '+str(chkwvn)+' cm^-1',flush=1)
    for date in dlist:
        raddfile='diag_'+sensor+'_'+loop+'.'+str(date)+'.nc4'
        infile1=archdir+'/'+str(date)+'/'+raddfile
    
        if (os.path.exists(infile1)):
            #print('Processing Radfile: %s' %(raddfile),flush=1)
            ds1=xa.open_dataset(infile1)
            npts=int(ds1.nobs.size/ds1.nchans.size)
            nchs=ds1.nchans.size
            #ds1=ds1.swap_dims({"nchans":"wavenumber"}) #replace the dimension of channel by channel indices
            #wavelength=1e+04/ds1.wavenumber
        else:
            #print('%s is not existing'%(raddfile),flush=1)
            continue
        
        # Observation lat/lon
        rlat1=np.reshape(ds1.Latitude.values,(npts,nchs))[:,0]
        rlon1=np.reshape(ds1.Longitude.values,(npts,nchs))[:,0]
        qcflags=np.reshape(ds1.QC_Flag.values,(npts,nchs))
        # obs1=np.reshape(ds1.Observation.values,(npts,nchs))
        sim1=np.reshape(ds1.Simulated_Tb.values,(npts,nchs))
        clr1=np.reshape(ds1.Clearsky_Tb.values,(npts,nchs))
        varinv1=np.reshape(ds1.Inverse_Observation_Error.values,(npts,nchs))
        # omb_bc1=np.reshape(ds1.Obs_Minus_Forecast_adjusted.values,(npts,nchs))
        omb_nbc1=np.reshape(ds1.Obs_Minus_Forecast_unadjusted.values,(npts,nchs))
        obs1=omb_nbc1+sim1
        tmpds=xa.Dataset({'rlon1':(['obsloc'],rlon1),
                          'rlat1':(['obsloc'],rlat1),
                          'qcflag':(['obsloc','wavenumber'],qcflags),
                          'tb_obs':(['obsloc','wavenumber'],obs1),
                          'tb_sim':(['obsloc','wavenumber'],sim1),
                          'tb_clr':(['obsloc','wavenumber'],clr1),
                          'varinv':(['obsloc','wavenumber'],varinv1),
                          'omb_nbc':(['obsloc','wavenumber'],omb_nbc1)},
                         coords={'obsloc':np.arange(npts),
                                 'wavenumber':ds1.wavenumber.values})
        tmpds=tmpds.sel(wavenumber=chkwvn)
        
        if (date==str(sdate)):
            ds_chk=tmpds
        else:
            ds_chk=xa.concat((ds_chk,tmpds),dim='obsloc')

    print('Finish loading data for channel with wavenumber '+str(chkwvn)+' cm^-1',flush=1)
    #ds_chk=ds_all#.sel(wavenumber=chkwvn)
    tb_sim=ds_chk.tb_sim
    tb_clr=ds_chk.tb_clr
    tb_obs=ds_chk.tb_obs
    varinv=ds_chk.varinv
    info_tmp=ds_sensorinfo.sel(wavenumber=chkwvn)
    nuchrad=info_tmp.nuchan.values
    iuserad=info_tmp.iuse.values
    obserr=info_tmp.obserr.values

    omb=ds_chk.omb_nbc
    aereff_fg=tb_sim-tb_clr
    aereff_obs=tb_obs-tb_clr
    aereff=0.5*abs(aereff_fg)+0.5*abs(aereff_obs)
    
    good_msk=(ds_chk.qcflag==0.)
    gross_msk=(ds_chk.qcflag==3.)
    cld_msk=(ds_chk.qcflag==7.)
    tzr_msk=(ds_chk.qcflag==10.)
    aer_msk=(ds_chk.qcflag==13.)
    sfcir_msk=(ds_chk.qcflag==53.)
    bust_msk=(ds_chk.qcflag==55.)
    #bust_msk=(aer_msk)&((abs(omb)>3)&(abs(omb)>1.8*aereff))

    ori_msk=((good_msk)|(aer_msk)|(gross_msk)|(sfcir_msk)|(tzr_msk))
    ori_total=np.count_nonzero(ori_msk)
    final_qc_msk=(good_msk)|(aer_msk)

    halfbin=0.5*binsize
    hist_x_edge=np.arange(-1*halfbin,50.+binsize,binsize)
    bin_center=(hist_x_edge+halfbin)[:-1]
    
    omb_mean2=np.zeros_like(bin_center,dtype='float')
    omb_sd2=np.zeros_like(bin_center,dtype='float')
    sd_noqc=omb[ori_msk==1].std()
    sd_qc=omb[final_qc_msk==1].std()    

    for i in np.arange(omb_mean2.size):
        lb_aereff=hist_x_edge[i]
        ub_aereff=hist_x_edge[i+1]
        tmpmsk2=(final_qc_msk)&((aereff>=lb_aereff)&(aereff<ub_aereff))
        omb_mean2[i]=omb[tmpmsk2==1].mean()
        omb_sd2[i]=omb[tmpmsk2==1].std()

    tmpsd2=xa.where(np.isnan(omb_sd2),-999,omb_sd2)
    gtmask=(tmpsd2>obserr)
    if (not any(gtmask)):
        Aeff_1_idx=0
    else:
        Aeff_1_idx=np.where(gtmask)[0].min()
    Aeff_2_idx=tmpsd2.argmax()
    aereff_1=bin_center[Aeff_1_idx]
    aereff_2=bin_center[Aeff_2_idx]
    sd_min=omb_sd2[Aeff_1_idx]
    sd_max=omb_sd2[Aeff_2_idx]
    if (sd_max > obserr):
        aero_sensitive=1
    else:
        aero_sensitive=0

    tmpdf=pd.DataFrame(np.array([[chkwvn,nuchrad,iuserad,obserr,aereff_1,aereff_2,sd_min,sd_max,sd_noqc,sd_qc,aero_sensitive]]),
                       columns=['wavenumber','nuchan','iuse','SD_o','Aeff_1','Aeff_2','SD_min','SD_max','SD_noqc','SD_qc','Aer_sen'])

    if (icount==0):
        df_all=tmpdf
    else:
        df_all=pd.concat((df_all,tmpdf))
    icount+=1

df_all.to_csv(savedir+'/'+sensor+'_'+str(nchs)+'_stats_new_seg'+str(seg)+'.'+version+'.csv')

#pr.disable()
#s = io.StringIO()
#sortby = SortKey.CUMULATIVE
#ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
#ps.print_stats()
#print(s.getvalue())

# df_all.to_excel(savedir+'/'+sensor+'_616_stats.xlsx')
# ds_all=ds_all.assign({'obserr':(['wavenumber'],err_array)})
# ds_all=ds_all.sel(wavenumber=spectral_range)