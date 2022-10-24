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
from plot_utils import setupax_2dmap, plt_x2y, set_size
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
spectral_range=slice(750,1300)
loop='ges' #ges,anl
useqc=-1
usemsk=2
plt_sp=1 # plot spectrum
fill_std=0

area='Glb'
minlon, maxlon, minlat, maxlat, crosszero, cyclic=setarea.setarea(area)
print(minlat,maxlat,minlon,maxlon,crosszero,cyclic)
if (area=='Glb'):
   minlon=-180. ; maxlon=180.
cornll=[minlat,maxlat,minlon,maxlon]

if (loop=='ges'):
   ylb='Total BC [K]'
else:
   ylb='Total BC [K]'

if (useqc==-2):
    qcflg='noqc'
elif (useqc==-1):
    qcflg='qc'
else:
    qcflg='qc%s'%(useqc)

rmflier=1
wateronly=0
if (wateronly):
    waterflg='water'
else:
    waterflg='all'
if (usemsk==-1):
    mskflg='omsk'
elif (usemsk==0):
    mskflg='msk0'
elif (usemsk==1):
    mskflg='msk1'
elif (usemsk==2):
    mskflg='imsk'
elif (usemsk==3):
    mskflg='iaer'


# Data path setup
lutpath=rootpath+'/AlbanyWork/Prospectus/Experiments/AeroObsStats/SD_LUT'
outpath=rootpath+'/AlbanyWork/Prospectus/Experiments/HazyDA/Images/DiagFiles/rad'
archdir0=rootarch+'/'+explist[0]
archdir1=rootarch+'/'+explist[1]

imgsavpath=outpath+'/spectrum/bias/'+area
if ( not os.path.exists(imgsavpath) ):
    os.makedirs(imgsavpath)

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
        print('Processing Radfile: %s' %(raddfile1),flush=1)
        if (dates_count==0):
           ds0=xa.open_dataset(infile0)
           ds0=ds0.swap_dims({"nchans":"wavenumber"})
           chkwvn_list=ds0.wavenumber.sel(wavenumber=spectral_range)[ds0.use_flag.sel(wavenumber=spectral_range)==1]
        dates_count+=1
    else:
        print('%s is not existing'%(raddfile1),flush=1)
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

cnts0=ds_all0.obsloc.size
cnts1=ds_all1.obsloc.size
ds_all0=ds_all0.assign_coords(obsloc=np.arange(cnts0))
ds_all1=ds_all1.assign_coords(obsloc=np.arange(cnts1))

wrkdata0=ds_all0.omb_nbc-ds_all0.omb_bc
wrkdata1=ds_all1.omb_nbc-ds_all1.omb_bc

mask0=~np.isnan(ds_all0.rlon)
mask1=~np.isnan(ds_all1.rlon)

if (useqc==-2):
    pass
elif (useqc==-1):
    mask0=(mask0)&((ds_all0.qcflag==0)|(ds_all0.qcflag==13))
    mask1=(mask1)&((ds_all1.qcflag==0)|(ds_all1.qcflag==13))
else:
    mask0=(mask0)&((ds_all0.qcflag==useqc))
    mask1=(mask1)&((ds_all1.qcflag==useqc))

if (usemsk==-1):
    pltmsk0=mask0
    pltmsk1=mask1
elif (usemsk==0):
    pltmsk0=mask0
    pltmsk1=mask0
elif (usemsk==1):
    pltmsk0=mask1
    pltmsk1=mask1
elif (usemsk==2):
    pltmsk0=mask0&mask1
    pltmsk1=mask0&mask1
elif (usemsk==3):
    pltmsk0=mask0&mask1&(ds_all1.qcflag==13)
    pltmsk1=mask0&mask1&(ds_all1.qcflag==13)

print('plot counts: %s, %s' %(np.count_nonzero(pltmsk0),np.count_nonzero(pltmsk1)))

wvn=ds_all0.wavenumber.data
wvl=1e4/wvn
wvnlb='Wavenumber [$cm^{-1}$]'
wvllb='Wavelength [µm]'
prop_dict={'color'     :['b','r'],
           'line_style':[' ',' '],
           'line_width':[1.5,1.5],
           'marker'    :['v','v'],
           'mark_size' :[6.,6.],
           'fillstyle' :['full','full'],
           'legend'    :expnlist,
           }
tistr=''

if (plt_sp):            
   pltda_x0=xa.where(pltmsk0,wrkdata0,np.nan)
   pltda_x1=xa.where(pltmsk1,wrkdata1,np.nan)

   mean0=pltda_x0.mean(dim='obsloc',skipna=1)
   mean1=pltda_x1.mean(dim='obsloc',skipna=1)

   stdv0=pltda_x0.std(dim='obsloc',skipna=1)
   stdv1=pltda_x1.std(dim='obsloc',skipna=1)

   fig,ax=plt.subplots()
   set_size(axe_w,axe_h,ax=ax,b=0.25)
   yaxlb='Mean %s' %(ylb)

   pltds=xa.Dataset({'pltdata':(['exps','channels'],xa.concat((mean0,mean1),dim='exps').data),
                     'pltstdv':(['exps','channels'],xa.concat((stdv0,stdv1),dim='exps').data),
                     },coords={'exps':explist,'channels':chkwvn_list.data})
   plt_x2y(pltds,yaxlb,wvn,wvnlb,wvl,wvllb,prop_dict,tistr,0,[],fill_std=fill_std,plot_diff=1,ax=ax)
   
   if (fsave):
      outname='%s/%s_%s_%s_%s.%s_%s.%s' %(imgsavpath,sensor,loop,qcflg,mskflg,sdate,edate,ffmt)
      print(outname,flush=1)
      fig.savefig(outname,dpi=quality)
      plt.close()

