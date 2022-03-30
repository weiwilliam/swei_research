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
axe_w=7 ; axe_h=3 ; quality=300

# Projection setting
proj=ccrs.PlateCarree(globe=None)

# Plotting setup
sdate=2020061000
edate=2020071018
aertag='Dust'
hint=6
exp='hazyda_aero'
expn='AER'
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
archpath=rootarch+'/'+exp
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
    infile1=archpath+'/'+str(date)+'/'+raddfile

    if (os.path.exists(infile1)):
        print('Processing Radfile: %s' %(raddfile))
        ds1=xa.open_dataset(infile1)
        npts=int(ds1.nobs.size/ds1.nchans.size)
        nchs=ds1.nchans.size
        ds1=ds1.assign_coords(nuchan=('wavenumber',ds1.wavenumber.data))
        ds1=ds1.swap_dims({"nchans":"wavenumber"}) #replace the dimension of channel by channel indices
        wavelength=1e+04/ds1.wavenumber
        chkwvn_list=ds1.wavenumber.sel(wavenumber=spectral_range)[ds1.use_flag.sel(wavenumber=spectral_range)==1]
        
        rlat1=np.reshape(ds1.Latitude.values,(npts,nchs))
        rlon1=np.reshape(ds1.Longitude.values,(npts,nchs))
        qcflags=np.reshape(ds1.QC_Flag.values,(npts,nchs))
        obs1=np.reshape(ds1.Observation.values,(npts,nchs))
        sim1=np.reshape(ds1.Simulated_Tb.values,(npts,nchs))
        clr1=np.reshape(ds1.Clearsky_Tb.values,(npts,nchs))
        tmpds=xa.Dataset({'rlon1':(['obsloc'],rlon1[:,0]),
                          'rlat1':(['obsloc'],rlat1[:,0]),
                          'qcflag':(['obsloc','wavenumber'],qcflags),
                          'tb_obs':(['obsloc','wavenumber'],obs1),
                          'tb_sim':(['obsloc','wavenumber'],sim1),
                          'tb_clr':(['obsloc','wavenumber'],clr1)},
                         coords={'obsloc':np.arange(npts),
                                 'wavenumber':ds1.wavenumber.values})
        
        tmpds=tmpds.sel(wavenumber=chkwvn_list)
        
        qc0cnts=np.count_nonzero((tmpds.qcflag==0),axis=0)
        qc13cnts=np.count_nonzero((tmpds.qcflag==13),axis=0)
        
    else:
        print('%s is not existing'%(raddfile))
        if ('tmpds' in locals()):
            qc0cnts=xa.zeros_like(tmpds.tb_obs[0,:])
            qc13cnts=xa.zeros_like(tmpds.tb_obs[0,:])
            qc0cnts[:]=np.nan
            qc13cnts[:]=np.nan
        else:
            print('postpone the starting date')
            continue
    
    qc_cnts=xa.Dataset({'qc0cnts':(['wavenumber'],qc0cnts.data),
                        'qc13cnts':(['wavenumber'],qc13cnts.data)},
                       coords={'wavenumber':chkwvn_list})

    # Observation lat/lon
    if (date==str(sdate)):
        qccnts_all=qc_cnts
    else:
        qccnts_all=xa.concat((qccnts_all,qc_cnts),dim='dates')

qccnts_all=qccnts_all.assign_coords(dates=('dates',dates))

savedir=outpath+'/'+expn+'/timeseries/'+aertag
if ( not os.path.exists(savedir) ):
    os.makedirs(savedir)
# for chkwvn in [906.25]:
for chkwvn in chkwvn_list:
    qc_chk=qccnts_all.sel(wavenumber=chkwvn)
    
    fig,ax=plt.subplots()
    set_size(axe_w,axe_h,b=0.2)
    ax.plot_date(dates,qc_chk.qc0cnts,'.-',color='tab:blue',linewidth=0.8)
    ax.plot_date(dates,qc_chk.qc13cnts,'.-',color='tab:red',linewidth=0.8)
    
    ax.xaxis.set_major_locator(loc)
    ax.xaxis.set_major_formatter(formatter)
    ax.xaxis.set_tick_params(rotation=20)
    
    tistr=('%s (%.2f $cm^{-1}$)' %(sensor,chkwvn))
    ax.set_title(tistr,loc='left')
    ax.set_xlabel('Dates')
    ax.set_ylabel('Counts')
            
    if (fsave):
        fname=('TS_%s_%s_%.2f.%s' %(area,sensor,chkwvn,ffmt))
        fig.savefig(savedir+'/'+fname,dpi=quality)
        plt.close()

savedir=outpath+'/'+expn+'/spec/'+aertag
if ( not os.path.exists(savedir) ):
    os.makedirs(savedir)

wvn=qccnts_all.wavenumber.values
wvl=1e4/wvn
wvnlb='Wavenumber [$cm^{-1}$]' 
wvllb='Wavelength [µm]'
colorlst=['black','tab:blue','tab:red']
lstylelst=[' ',' ',' ']
mrklst=['o','^','x']
lglst=['Total','Passed','Aer-affected']
yaxlb='Number of observations'

datessum=qccnts_all.sum(dim='dates')
total=datessum.qc0cnts+datessum.qc13cnts
qc0sum=datessum.qc0cnts
qc13sum=datessum.qc13cnts
   
pltda=xa.concat((total,qc0sum),dim='lines')
pltda=xa.concat((pltda,qc13sum),dim='lines').T
tistr=''
        
fig,ax=plt.subplots()
set_size(axe_w,axe_h,ax=ax,b=0.25)
plt_x2y(pltda,yaxlb,wvn,wvnlb,wvl,wvllb,colorlst,lstylelst,mrklst,tistr,lglst,0,[],ax=ax)

fname=('Spec_%s_%s.%s-%s.%s' %(area,sensor,spectral_range.start,spectral_range.stop,ffmt))
if (fsave):
    fig.savefig(savedir+'/'+fname,dpi=quality)
    plt.close()
