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
elif (os_name=='Linux'):
    rootpath='/scratch2/BMC/gsd-fv3-dev/Shih-wei.Wei'
    rootarch='/scratch2/BMC/gsd-fv3-dev/Shih-wei.Wei/ResearchData'
    rootgit='/home/Shih-wei.Wei/research'

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
mpl.rc('legend',fontsize='large')
fsave=1 ; ffmt='png' ; ptsize=4
axe_w=3 ; axe_h=3 ; quality=300

# Projection setting
proj=ccrs.PlateCarree(globe=None)

# Plotting setup
sdate=2020061000
edate=2020092118
aertype='All'
hint=6
exp='AerObserver'
sensor='iasi_metop-a'
spectral_range=slice(700,1300)
loop='ges' #ges,anl
binsize=0.1

# Data path setup
archpath=rootarch+'/Prospectus/AeroObsStats/nc_diag'
fixpath=rootgit+'/GSI_exps/fix'
outpath=rootpath+'/AlbanyWork/Prospectus/Experiments/AeroObsStats'
archdir=archpath+'/'+exp
print(archpath)
print(fixpath)
print(outpath)
print(archdir)
savedir=outpath+'/'+exp+'_newQC/SD_LUT/'+aertype
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
    print('Processing the first Radfile: %s' %(raddfile))
    ds1=xa.open_dataset(infile1)
    npts=int(ds1.nobs.size/ds1.nchans.size)
    nchs=ds1.nchans.size
    ds1=ds1.assign_coords(nuchan=('wavenumber',ds1.wavenumber))
    ds1=ds1.swap_dims({"nchans":"wavenumber"}) #replace the dimension of channel by channel indices
    wavelength=1e+04/ds1.wavenumber
    chkwvn_list=ds1.wavenumber.values
    ds_sensorinfo=xa.Dataset({'obserr':(['wavenumber'],err_array),
                              'nuchan':(['wavenumber'],nuch_array),
                              'iuse':(['wavenumber'],iuse_array)},
                             coords={'wavenumber':ds1.wavenumber.values})
else:
    print('%s is not existing'%(raddfile))

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
        # chkwvn_list=ds1.wavenumber.sel(wavenumber=spectral_range)[ds1.use_flag.sel(wavenumber=spectral_range)==1]
        #usedchidx=np.where(ds1.iuse_rad==1)[0]
        #unusedchidx=np.where(ds1.iuse_rad==-1)[0]
        #wvldiff=abs(np.subtract(wavelength[usedchidx],chkwvl))
        #chkwvlidx=usedchidx[np.where(wvldiff==wvldiff.min())[0][0]]
        #print('Check wavelength: %.2f' %(wavelength[chkwvlidx]))
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
    # tmpds=tmpds.sel(wavenumber=chkwvn)
    
    if (date==str(sdate)):
        ds_all=tmpds
    else:
        ds_all=xa.concat((ds_all,tmpds),dim='obsloc')
        
icount=0
for chkwvn in chkwvn_list:
# for chkwvn in [962.5]:
    print('Processing for channel with wavenumber '+str(chkwvn)+' cm^-1')
    ds_chk=ds_all.sel(wavenumber=chkwvn)
    tb_sim=ds_chk.tb_sim
    tb_clr=ds_chk.tb_clr
    tb_obs=ds_chk.tb_obs
    varinv=ds_chk.varinv
    info_tmp=ds_sensorinfo.sel(wavenumber=chkwvn)
    nuchrad=info_tmp.nuchan.values
    iuserad=info_tmp.iuse.values
    obserr=np.sqrt(info_tmp.obserr.values)

    omb=tb_obs-tb_sim
    aereff_fg=tb_sim-tb_clr
    aereff_obs=tb_obs-tb_clr
    aereff=0.5*abs(aereff_fg)+0.5*abs(aereff_obs)
    
    qc0_msk=(ds_chk.qcflag==0.)
    # qc7_msk=(ds_chk.qcflag==7.)
    qc13_msk=(ds_chk.qcflag==13.)

    ori_msk=((qc0_msk)|(qc13_msk))
    ori_total=np.count_nonzero(ori_msk)
    # lowbt_qc=(used_msk)&(abs(omb)<30.)
    final_qc_msk=(ori_msk)&(~((abs(omb)>3)&(abs(omb)>1.8*aereff)))&(abs(omb)<30.)

    halfbin=0.5*binsize
    hist_x_edge=np.arange(-1*halfbin,50.+binsize,binsize)
    bin_center=(hist_x_edge+halfbin)[:-1]
    
#    omb_mean1=np.zeros_like(bin_center,dtype='float')
#    omb_sd1=np.zeros_like(bin_center,dtype='float')
#    counts1=np.zeros_like(bin_center,dtype='int')
    omb_mean2=np.zeros_like(bin_center,dtype='float')
    omb_sd2=np.zeros_like(bin_center,dtype='float')
#    counts2=np.zeros_like(bin_center,dtype='int')
    sd_noqc=omb[ori_msk==1].std()
    sd_qc=omb[final_qc_msk==1].std()    

    for i in np.arange(omb_mean2.size):
        lb_aereff=hist_x_edge[i]
        ub_aereff=hist_x_edge[i+1]
#        tmpmsk1=(ori_msk)&((aereff>=lb_aereff)&(aereff<ub_aereff))
#        omb_mean1[i]=omb[tmpmsk1==1].mean()
#        omb_sd1[i]=omb[tmpmsk1==1].std()
#        counts1[i]=np.count_nonzero(tmpmsk1)
        tmpmsk2=(final_qc_msk)&((aereff>=lb_aereff)&(aereff<ub_aereff))
        omb_mean2[i]=omb[tmpmsk2==1].mean()
        omb_sd2[i]=omb[tmpmsk2==1].std()
#        counts2[i]=np.count_nonzero(tmpmsk2)

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

    tmpdf=pd.DataFrame(np.array([[chkwvn,nuchrad,iuserad,obserr,aereff_1,aereff_2,sd_min,sd_max,sd_noqc,sd_qc]]),
                       columns=['wavenumber','nuchan','iuse','SD_o','Aeff_1','Aeff_2','SD_min','SD_max','SD_noqc','SD_qc'])

    if (icount==0):
        df_all=tmpdf
    else:
        df_all=pd.concat((df_all,tmpdf))
    icount+=1

df_all.to_csv(savedir+'/'+sensor+'_'+str(nchs)+'_stats.csv')
# df_all.to_excel(savedir+'/'+sensor+'_616_stats.xlsx')

# ds_all=ds_all.assign({'obserr':(['wavenumber'],err_array)})
# ds_all=ds_all.sel(wavenumber=spectral_range)
