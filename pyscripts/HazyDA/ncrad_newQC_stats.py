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
machine='S4'
if (machine=='MBP'):
    rootpath='/Users/weiwilliam'
    rootarch='/Volumes/WD2TB/ResearchData'
elif (machine=='Desktop'):
    rootpath='F:\GoogleDrive_NCU\Albany'
    rootarch='F:\ResearchData'
    rootgit='F:\GitHub\swei_research'
elif (machine=='S4'):
    rootarch='/data/users/swei/ResearchData/Prospectus/AeroObsStats/nc_diag'
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
exp='aero_v2qc'
sensor='iasi_metop-a'
spectral_range=slice(700,1300)
loop='ges' #ges,anl
binsize=0.1
version='v5'

# Data path setup
#archpath=rootarch+'/AeroObsStats/OUTPUT'
fixpath=rootgit+'/GSI_exps/fix'
outpath=rootpath+'/AlbanyWork/Prospectus/Experiments/AeroObsStats'
archdir=rootarch+'/'+exp
#print(archpath)
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
    ds_sensorinfo=xa.Dataset({'obserr':(['wavenumber'],err_array),
                              'nuchan':(['wavenumber'],nuch_array),
                              'iuse':(['wavenumber'],iuse_array),
                              'ich':(['wavenumber'],np.arange(nchs)+1),
                              },
                             coords={'wavenumber':ds1.wavenumber.values})
else:
    print('%s is not existing'%(raddfile),flush=1)

#print('Processing for channel with wavenumber '+str(chkwvn)+' cm^-1',flush=1)
for date in dlist:
    raddfile='diag_'+sensor+'_'+loop+'.'+str(date)+'.nc4'
    infile1=archdir+'/'+str(date)+'/'+raddfile

    if (os.path.exists(infile1)):
        print('Processing Radfile: %s' %(raddfile),flush=1)
        ds1=xa.open_dataset(infile1)
        npts=int(ds1.nobs.size/ds1.nchans.size)
        nchs=ds1.nchans.size
        chkwvn_list=ds1.wavenumber[ds1.use_flag==1].data
        #ds1=ds1.swap_dims({"nchans":"wavenumber"}) #replace the dimension of channel by channel indices
        #wavelength=1e+04/ds1.wavenumber
    else:
        print('%s is not existing'%(raddfile),flush=1)
        continue
    
    # Observation lat/lon
    rlat1=np.reshape(ds1.Latitude.values,(npts,nchs))[:,0]
    rlon1=np.reshape(ds1.Longitude.values,(npts,nchs))[:,0]
    qcflags=np.reshape(ds1.QC_Flag.values,(npts,nchs))
    # obs1=np.reshape(ds1.Observation.values,(npts,nchs))
    sim1=np.reshape(ds1.Simulated_Tb.values,(npts,nchs))
    clr1=np.reshape(ds1.Clearsky_Tb.values,(npts,nchs))
    #varinv1=np.reshape(ds1.Inverse_Observation_Error.values,(npts,nchs))
    # omb_bc1=np.reshape(ds1.Obs_Minus_Forecast_adjusted.values,(npts,nchs))
    omb_nbc1=np.reshape(ds1.Obs_Minus_Forecast_unadjusted.values,(npts,nchs))
    obs1=omb_nbc1+sim1
    aereff_fg=sim1-clr1
    aereff_obs=obs1-clr1
    aereff=0.5*abs(aereff_fg)+0.5*abs(aereff_obs)
    tmpds=xa.Dataset({'rlon1':(['obsloc'],rlon1),
                      'rlat1':(['obsloc'],rlat1),
                      'qcflag':(['obsloc','wavenumber'],qcflags),
                      'SD_o':(['wavenumber'],ds_sensorinfo.obserr.data),
                      'aereff':(['obsloc','wavenumber'],aereff),
                      'omb_nbc':(['obsloc','wavenumber'],omb_nbc1)},
                     coords={'obsloc':np.arange(npts),
                             'wavenumber':ds1.wavenumber.data})
    tmpds=tmpds.sel(wavenumber=chkwvn_list)
    
    if (date==str(sdate)):
        ds_chk=tmpds
    else:
        ds_chk=xa.concat((ds_chk,tmpds),dim='obsloc')

halfbin=0.5*binsize
hist_x_edge=np.arange(-1*halfbin,50.+binsize,binsize)
bin_center=(hist_x_edge+halfbin)[:-1]

#omb=ds_chk.omb_nbc

#good_msk=(ds_chk.qcflag==0.)
#gross_msk=(ds_chk.qcflag==3.)
#cld_msk=(ds_chk.qcflag==7.)
#tzr_msk=(ds_chk.qcflag==10.)
#aer_msk=(ds_chk.qcflag==13.)
#sfcir_msk=(ds_chk.qcflag==53.)
#bust_msk=(ds_chk.qcflag==55.)
#aercld_msk=(ds_chk.qcflag==57.)
##bust_msk=(aer_msk)&((abs(omb)>3)&(abs(omb)>1.8*aereff))

##ori_msk=((good_msk)|(aer_msk)|(gross_msk)|(sfcir_msk)|(tzr_msk)|(aercld_msk)|(bust_msk))
#ori_msk=(~cld_msk)
#ori_total=np.count_nonzero(ori_msk)
#final_qc_msk=(good_msk)|(aer_msk)

print('whole period of data is loaded', flush=1)

ds_sinfo_chk=ds_sensorinfo.sel(wavenumber=chkwvn_list)
df0=ds_chk.to_dataframe()
df0=df0.reset_index()
preqc_filter=((df0['qcflag']!=7.0))
preqc_df=df0.loc[preqc_filter,:]
aftqc_filter=((df0['qcflag']==0.0)|(df0['qcflag']==13.0))
aftqc_df=df0.loc[aftqc_filter,:]

print('preqc_df and aftqc_df have been created',flush=1)

aftqc_df['Ae_bin']=pd.cut(aftqc_df['aereff'],bins=hist_x_edge,labels=bin_center)

preqc_grp=preqc_df.groupby(['wavenumber']).agg({'omb_nbc':['std']}).reset_index()
aftqc_grp=aftqc_df.groupby(['wavenumber']).agg({'omb_nbc':['std']}).reset_index()
tmp_grp=aftqc_df.groupby(['wavenumber','Ae_bin']).agg({'omb_nbc':['std'],'SD_o':['min']}).reset_index()

preqc_grp.columns=preqc_grp.columns.droplevel(1)
preqc_grp=preqc_grp.rename(columns={'omb_nbc':'SD_noqc'})
aftqc_grp.columns=aftqc_grp.columns.droplevel(1)
aftqc_grp=aftqc_grp.rename(columns={'omb_nbc':'SD_qc'})
tmp_grp.columns=tmp_grp.columns.droplevel(1)
tmp_grp=tmp_grp.rename(columns={'omb_nbc':'SD_qc'})

tmpdf=tmp_grp.loc[tmp_grp.SD_qc>tmp_grp.SD_o,:]
ae1idx=tmpdf.groupby(['wavenumber'])['Ae_bin'].transform(min)==tmpdf['Ae_bin']
Aeff1_df=tmpdf[ae1idx].set_index('wavenumber')
Aeff1_df=Aeff1_df.rename(columns={'Ae_bin':'Aeff_1','SD_qc':'SD_min'})

ae2idx=tmp_grp.groupby(['wavenumber'])['SD_qc'].transform(max)==tmp_grp['SD_qc']
Aeff2_df=tmp_grp[ae2idx].set_index('wavenumber')
Aeff2_df=Aeff2_df.rename(columns={'Ae_bin':'Aeff_2','SD_qc':'SD_max'})
print('Aeff1_df and Aeff2_df have been created',flush=1)

output_df=pd.concat((Aeff1_df,Aeff2_df),axis=1).sort_index()
output_df=output_df[['Aeff_1','Aeff_2','SD_min','SD_max']]
output_df=output_df.reset_index()
output_df['nuchan']=ds_sinfo_chk.nuchan.data
output_df['iuse']  =ds_sinfo_chk.iuse.data
output_df['ich']   =ds_sinfo_chk.ich.data
output_df['SD_o']  =ds_sinfo_chk.obserr.data
output_df['SD_noqc']=preqc_grp['SD_noqc']
output_df['SD_qc']=aftqc_grp['SD_qc']
output_df['Aer_sen']=(output_df['SD_max']>output_df['SD_o']).astype(int)
output_df=output_df[['wavenumber','ich','nuchan','iuse','SD_o','Aeff_1','Aeff_2','SD_min','SD_max','SD_noqc','SD_qc','Aer_sen']]
print('output_df has been created',flush=1)

#
outfile=savedir+'/'+sensor+'_'+str(nchs)+'_stats.'+version+'.csv'
output_df.to_csv(outfile)
print(outfile+' has been created',flush=1)
