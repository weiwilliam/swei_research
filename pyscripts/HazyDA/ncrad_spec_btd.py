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
from plot_utils import setupax_2dmap, plt_1exp_x2y, set_size
from utils import ndate,setup_cmap
from gsi_ncdiag import read_rad_ncdiag

tlsize=12 ; lbsize=10
mpl.rc('axes', titlesize=tlsize,labelsize=lbsize)
mpl.rc('xtick',labelsize=lbsize)
mpl.rc('ytick',labelsize=lbsize)
mpl.rc('legend',fontsize='small')
fsave=1 ; ffmt='png' ; ptsize=4
axe_w=8 ; axe_h=3 ; quality=300

# Plotting setup
sdate=2020061000
edate=2020071018
aertype='Dust'
hint=6
exp='hazyda_aero'
expn=['AER']
sensor='iasi_metop-a'
loop='ges' #ges,anl
plt_sp=1 # plot spectrum
fill_std=0

area='Glb'
minlon, maxlon, minlat, maxlat, crosszero, cyclic=setarea.setarea(area)
print(minlat,maxlat,minlon,maxlon,crosszero,cyclic)
if (area=='Glb'):
   minlon=-180. ; maxlon=180.
cornll=[minlat,maxlat,minlon,maxlon]

# Data path setup
lutpath=rootpath+'/AlbanyWork/Prospectus/Experiments/AeroObsStats/SD_LUT'
outpath=rootpath+'/AlbanyWork/Prospectus/Experiments/HazyDA/Images/DiagFiles/rad'
archdir=rootarch+'/'+exp

imgsavpath=outpath+'/spectrum/btd/'+area
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
    raddfile='diag_'+sensor+'_'+loop+'.'+str(date)+'.nc4'
    infile=archdir+'/'+str(date)+'/'+raddfile

    if (os.path.exists(infile)):
        print('Processing Radfile: %s' %(raddfile),flush=1)
        if dates_count==0:
           ds=xa.open_dataset(infile)
           ds=ds.swap_dims({"nchans":"wavenumber"})
           chkwvn_list=ds.wavenumber[ds.use_flag==1]
        dates_count+=1
    else:
        print('%s is not existing'%(raddfile),flush=1)
        continue

    tmpds=read_rad_ncdiag(infile,chkwvn=chkwvn_list)
      
    if (date==str(sdate)):
        ds_all=tmpds
    else:
        ds_all=xa.concat((ds_all,tmpds),dim='obsloc')

cnts=ds_all.obsloc.size
ds_all=ds_all.assign_coords(obsloc=np.arange(cnts))
btd=ds_all.tb_sim-ds_all.tb_clr
pltmsk=(ds_all.qcflag==13.)|(ds_all.qcflag==57.)

#print('plot counts: %s' %(np.count_nonzero(pltmsk)))

wvn=ds_all.wavenumber.data
wvl=1e4/wvn
wvnlb='Wavenumber [$cm^{-1}$]'
wvllb='Wavelength [µm]'
prop_dict={'color'     :['g','r'],
           'line_style':['-','-'],
           'line_width':[1.5,1.5],
           'marker'    :['o','v'],
           'mark_size' :[3.,3.],
           'fillstyle' :['none','none'],
           'legend'    :expn,
           }
tistr=''

if (plt_sp):            
   pltda_x=xa.where(pltmsk,btd,np.nan)
   mean=pltda_x.mean(dim='obsloc',skipna=1)
   stdv=pltda_x.std(dim='obsloc',skipna=1)

   fig,ax=plt.subplots()
   set_size(axe_w,axe_h,ax=ax,b=0.25,r=0.9)
   yaxlb='%s' %('Mean BTD [K]')

   pltds=xa.Dataset({'pltdata':(['channels'],mean.data),
                     'pltstdv':(['channels'],stdv.data),
                     },coords={'channels':wvn})
   plt_1exp_x2y(pltds,yaxlb,wvn,wvnlb,wvl,wvllb,prop_dict,tistr,0,[],fill_std=fill_std,ax=ax)
   
   if (fsave):
      outname='%s/%s_%s_%s_btd.%s_%s.%s' %(imgsavpath,exp,sensor,loop,sdate,edate,ffmt)
      print(outname,flush=1)
      fig.savefig(outname,dpi=quality)
      plt.close()

