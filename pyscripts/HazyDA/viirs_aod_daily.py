#!/usr/bin/env python3
import os,sys, platform
machine='S4'
if (machine=='MBP'):
    rootpath='/Users/weiwilliam'
    rootarch='/Volumes/WD2TB/ResearchData'
elif (machine=='Desktop'):
    rootpath='F:\GoogleDrive_NCU\Albany'
    rootarch='F:\ResearchData'
    rootgit='F:\GitHub\swei_research'
elif (machine=='S4'):
    rootpath='/data/users/swei/AlbanyWork/Prospectus/Experiments/HazyDA/Images'
    rootarch='/data/users/swei'
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
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mpcrs
import setuparea as setarea
from plot_utils import setupax_2dmap, plt_x2y, set_size
from utils import ndate,setup_cmap
from datetime import datetime, timedelta
import cartopy.crs as ccrs

tlsize=12 ; txsize=12
mpl.rc('axes', titlesize=tlsize,labelsize=txsize)
mpl.rc('xtick',labelsize=txsize)
mpl.rc('ytick',labelsize=txsize)
mpl.rc('legend',fontsize='small')
fsave=1 ; ffmt='png' ; ptsize=4
axe_w=6 ; axe_h=3 ; quality=300
tkfreq=2
minussign=u'\u2212'

# Projection setting
proj=ccrs.PlateCarree(globe=None)

viirs_aod_path='/data/users/swei/Dataset/VIIRS_NPP/AOD_daily_1by1'
viirs_1x1_file='AERDB_D3_VIIRS_SNPP.A2020192.011.2021099192750.nc'

sdate='2020061000'
edate='2020071000'
hint=24
var='Aerosol_Optical_Thickness_550_Land_Ocean_Mean'
area='Glb'

minlon, maxlon, minlat, maxlat, crosszero, cyclic=setarea.setarea(area)
print(area,minlat,maxlat,minlon,maxlon,crosszero)
if (area=='Glb'):
   minlon=-180. ; maxlon=180.
else:
   minlon=(minlon+180)%360-180
   maxlon=(maxlon+180)%360-180
cornerll=[minlat,maxlat,minlon,maxlon]

outputpath=rootpath
imgsavpath=outputpath+'/VIIRS/L3AOD_2d/'+area
if ( not os.path.exists(imgsavpath) ):
   os.makedirs(imgsavpath)

cbori='horizontal' #vertical, horizontal
if (cbori=='vertical'):
   cb_frac=0.025
   cb_pad=0.06
elif (cbori=='horizontal'):
   cb_frac=0.04
   cb_pad=0.1

date1 = pd.to_datetime(sdate,format='%Y%m%d%H')
date2 = pd.to_datetime(edate,format='%Y%m%d%H')
delta = timedelta(hours=24)
dates = pd.date_range(start=date1, end=date2, freq=delta)

i=0
for file in sorted(os.listdir(viirs_aod_path)):
   tmp_data_time=pd.to_datetime(file.split('.')[1][1:],format='%Y%j')
   if tmp_data_time in dates:
      i+=1
      print(file,flush=1)
      ds=xa.open_dataset(os.path.join(viirs_aod_path,file))
   else:
      continue
   tmp=ds[var]
   if i==1:
      aod=tmp
   else:
      aod=xa.concat((aod,tmp),dim='time')

cnlvs=np.array((0., 0.05, 0.1, 0.15, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 2.5))
clridx=np.array((0,2,4,7,8,9,10,12,14,15,16,17,18))
clrmap=setup_cmap('precip3_16lev',clridx)
aer_norm = mpcrs.BoundaryNorm(cnlvs,len(clridx),extend='both')
cblb='VIIRS L3 AOD at 550 nm'

fig,ax,gl=setupax_2dmap(cornerll,area,proj,lbsize=txsize)
set_size(axe_w,axe_h,b=0.15,l=0.05,r=0.95,t=0.95)
pltdata=aod.mean(dim='time',skipna=1)
cn=ax.contourf(pltdata.Longitude_1D,pltdata.Latitude_1D,pltdata,levels=cnlvs,
               cmap=clrmap,norm=aer_norm,extend='both')
plt.colorbar(cn,orientation=cbori,fraction=cb_frac,
             pad=cb_pad,ticks=cnlvs[::tkfreq],aspect=40,label=cblb)
outname='%s/VIIRS_L3_AOD550.1x1.%s_%s.%s' %(imgsavpath,sdate,edate,ffmt)
if (fsave):
   print(outname,flush=1)
   fig.savefig(outname,dpi=quality)
   plt.close()



