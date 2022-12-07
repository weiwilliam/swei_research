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
import seaborn as sb
from datetime import datetime, timedelta
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mpcrs
import cartopy.crs as ccrs
import setuparea as setarea
from plot_utils import setupax_2dmap, set_size
from utils import ndate,setup_cmap
from gsi_ncdiag import read_rad_ncdiag

tlsize=12 ; lbsize=10
mpl.rc('axes', titlesize=tlsize,labelsize=lbsize)
mpl.rc('xtick',labelsize=lbsize)
mpl.rc('ytick',labelsize=lbsize)
mpl.rc('legend',fontsize='small')
fsave=1 ; ffmt='png' ; ptsize=4
axe_w=6 ; axe_h=3 ; quality=300

# Plotting setup
sdate=2020061000
edate=2020071018
aertype='Dust'
hint=6
explist=['hazyda_ctrl','hazyda_aero']
expnlist=['CTL','AER']
sensor='iasi_metop-a'
spectral_range=slice(700,1300)
loop='ges' #ges,anl
usebc=0
pltbx=1 # plot 2d histogram

area='Glb'
minlon, maxlon, minlat, maxlat, crosszero, cyclic=setarea.setarea(area)
print(minlat,maxlat,minlon,maxlon,crosszero,cyclic)
if (area=='Glb'):
   minlon=-180. ; maxlon=180.
cornll=[minlat,maxlat,minlon,maxlon]

if (usebc):
    bcflg='bc'
else:
    bcflg='nobc'

if (loop=='ges'):
   ylb='OMB [K]'
else:
   ylb='OMA [K]'
xlb='Wavenumber [$cm^{-1}$]'

# Data path setup
lutpath=rootpath+'/AlbanyWork/Prospectus/Experiments/AeroObsStats/SD_LUT'
outpath=rootpath+'/AlbanyWork/Prospectus/Experiments/HazyDA/Images/DiagFiles/rad'
archdir0=rootarch+'/'+explist[0]
archdir1=rootarch+'/'+explist[1]

savedir=outpath+'/boxplot/'+area
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
        os.path.exists(infile1)):
        print('Processing Radfile: %s' %(raddfile1))
        if (dates_count==0):
           ds0=xa.open_dataset(infile0)
           ds0=ds0.swap_dims({"nchans":"wavenumber"})
           chkwvn_list=ds0.wavenumber.sel(wavenumber=spectral_range)[ds0.use_flag.sel(wavenumber=spectral_range)==1]
        dates_count+=1
    else:
        print('%s is not existing'%(raddfile1))
        continue

    # Observation lat/lon from exp 0 (baseline)
    tmpds0=read_rad_ncdiag(infile0,chkwvn=chkwvn_list)
    tmpds1=read_rad_ncdiag(infile1,chkwvn=chkwvn_list)
      
    if (date==str(sdate)):
        ds_all0=tmpds0
        ds_all1=tmpds1
    else:
        ds_all0=xa.concat((ds_all0,tmpds0),dim='obsloc')
        ds_all1=xa.concat((ds_all1,tmpds1),dim='obsloc')

total_obscounts=ds_all1.obsloc.size
ds_all1=ds_all1.assign_coords(obsloc=np.arange(total_obscounts))

if (usebc):
    omb0=ds_all0.omb_bc
    omb1=ds_all1.omb_bc
else:
    omb0=ds_all0.omb_nbc
    omb1=ds_all1.omb_nbc

good_msk0  =(ds_all0.qcflag==0.)
gross_msk0 =(ds_all0.qcflag==3.)
cld_msk0   =(ds_all0.qcflag==7.)
tzr_msk0   =(ds_all0.qcflag==10.)
sfcir_msk0 =(ds_all0.qcflag==53.)
passed_msk0=(good_msk0)
# ori_msk=(good_msk)|(aer_msk)|(bust_msk)|(tzr_msk)

good_msk1  =(ds_all1.qcflag==0.)
gross_msk1 =(ds_all1.qcflag==3.)
cld_msk1   =(ds_all1.qcflag==7.)
tzr_msk1   =(ds_all1.qcflag==10.)
aer_msk1   =(ds_all1.qcflag==13.)
sfcir_msk1 =(ds_all1.qcflag==53.)
bust_msk1  =(ds_all1.qcflag==55.)
passed_msk1=(good_msk1)|(aer_msk1)
# ori_msk=(good_msk)|(aer_msk)|(bust_msk)|(tzr_msk)

pltmsk0=passed_msk0
pltmsk1=passed_msk1

if (pltbx):            
    pltda_x0=xa.where(pltmsk0,omb0,np.nan)
    pltda_x0=pltda_x0.assign_coords(exp=expnlist[0])
    pltda_x1=xa.where(pltmsk1,omb1,np.nan)
    pltda_x1=pltda_x1.assign_coords(exp=expnlist[1])

    chnlsize=pltda_x0.wavenumber.size
    if ( chnlsize > 30 ):
       chxint=6
       chnlseg=5
       chidxlist=np.ceil(np.linspace(0,chnlsize,chnlseg+1))

       for s in np.arange(chnlseg):
           sidx=int(chidxlist[s]); eidx=int(chidxlist[s+1])
           tmpda_x0=pltda_x0.sel(wavenumber=pltda_x0.wavenumber.data[sidx:eidx])
           tmpda_x1=pltda_x1.sel(wavenumber=pltda_x1.wavenumber.data[sidx:eidx])

           pltdf_x0=tmpda_x0.to_dataframe(ylb)
           pltdf_x1=tmpda_x1.to_dataframe(ylb)
    
           boxda=pd.concat((pltdf_x0,pltdf_x1))
           boxda=boxda.reset_index(level='wavenumber')
           boxda=boxda.rename(columns={'wavenumber':xlb})
           
           fig,ax=plt.subplots()
           set_size(axe_w,axe_h,b=0.18)
           ax=sb.boxplot(x=xlb,y=ylb,data=boxda,hue='exp',showfliers=False,whis=[5,95])
           ax.xaxis.set_ticks(np.arange(0,tmpda_x0.wavenumber.size,chxint))
           ax.xaxis.set_ticklabels(tmpda_x0.wavenumber.data[::chxint])
           ax.hlines(0,0,1,transform=ax.get_yaxis_transform(),colors='grey',linewidth=0.4)
           
           if (fsave):
               fname=('%s/BOX_%s_%s_%s.seg%s.%s'
                       %(savedir,sensor,loop,bcflg,s,ffmt))
               print(fname)
               fig.savefig(fname,dpi=quality)
               plt.close()
    else:
       chxint=5
       pltdf_x0=pltda_x0.to_dataframe(ylb)
       pltdf_x1=pltda_x1.to_dataframe(ylb)

       boxda=pd.concat((pltdf_x0,pltdf_x1))
       boxda=boxda.reset_index(level='wavenumber')
       boxda=boxda.rename({'wavenumber':xlb})

       fig,ax=plt.subplots()
       set_size(axe_w,axe_h,b=0.18)
       ax=sb.boxplot(x=xlb,y=ylb,data=boxda,hue='exp',showfliers=False,whis=[5,95])
       ax.xaxis.set_ticks(np.arange(0,pltda_x0.wavenumber.size,chxint))
       ax.xaxis.set_ticklabels(pltda_x0.wavenumber.data[::chxint])
       ax.hlines(0,0,1,transform=ax.get_yaxis_transform(),colors='grey',linewidth=0.4)

       if (fsave):
           fname=('%s/BOX_%s_%s_%s.allchnl.%s'
                   %(savedir,sensor,loop,bcflg,ffmt))
           print(fname)
           fig.savefig(fname,dpi=quality)
           plt.close()

