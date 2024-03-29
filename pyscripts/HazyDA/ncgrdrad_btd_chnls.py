#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 21:59:35 2019

@author: weiwilliam
"""
import sys, os, platform
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mpcrs
import xarray as xa
import pandas as pd

import setuparea as setarea
from plot_utils import setupax_2dmap, plt_x2y, set_size
from utils import setup_cmap,find_cnlvs
import cartopy.crs as ccrs

tlsize=10 ; txsize=10
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

sensorlist=['airs_aqua','amsua_aqua','amsua_metop-a','amsua_n15','amsua_n18',
            'amsua_n19','atms_npp','avhrr_metop-a','avhrr_n18','cris_npp','gmi_gpm',
            'hirs4_metop-a','hirs4_metop-b','hirs4_n19','iasi_metop-a','iasi_metop-b',
            'mhs_metop-a','mhs_metop-b','mhs_n18','mhs_n19','saphir_meghat',
            'seviri_m08','seviri_m10','sndrd1_g15','sndrd2_g15','sndrd3_g15',
            'sndrd4_g15','ssmis_f17','ssmis_f18']
#hsensorlist=['airs_aqua','iasi_metop-a','iasi_metop-b','cris_npp']
#lsensor1list=['hirs4_metop-a','hirs4_metop-b','hirs4_n19']
#lsensor2list=['sndrd1_g15','sndrd2_g15','sndrd3_g15','sndrd4_g15']
#lsensor3list=['avhrr_metop-a','avhrr_n18','seviri_m08','seviri_m10']

pltvar='btd_max'
varlb='BTD Maximum'
data_sdate=2020060106
data_edate=2020071018
check_sdate=2020060106
check_edate=2020071018
hint=6

sensor='iasi_metop-a'
chkwvn=962.5
degres=2.5
units='K'

exp='hazyda_aero'
expn='AER'

if check_sdate==data_sdate and check_edate==data_edate:
   grdtype='mean'
else:
   grdtype='time'

area='Glb'
minlon, maxlon, minlat, maxlat, crosszero, cyclic=setarea.setarea(area)
print(minlat,maxlat,minlon,maxlon,crosszero,cyclic)
if (area=='Glb'):
   minlon=-180. ; maxlon=180.
else:
   minlon=(minlon+180)%360-180
   maxlon=(maxlon+180)%360-180
cornerll=[minlat,maxlat,minlon,maxlon]

loop='ges' #ges,anl
if loop=='anl':
    tlstr='OMA'
elif loop=='ges':
    tlstr='OMF'

useqc=1
if (useqc):
   qcflg='qc'
else:
   qcflg='noqc'

cbori='horizontal' #vertical, horizontal
if (cbori=='vertical'):
   cb_frac=0.025
   cb_pad=0.06
elif (cbori=='horizontal'):
   cb_frac=0.04
   cb_pad=0.1

inpath=rootarch+'/archive/HazyDA/gridded_diag'
outputpath=rootpath+'/DiagFiles/gridded_rad'
imgsavpath=outputpath+'/spectrum/btd/'+area
if ( not os.path.exists(imgsavpath) ):
   os.makedirs(imgsavpath)

grdfile0='%s/%s_%s_%s_%s_btd_%.1fx%.1f.%s.%s_%s.nc' %(inpath,expn,sensor,loop,qcflg,degres,degres,grdtype,data_sdate,data_edate)

ds0=xa.open_dataset(grdfile0)
if grdtype == 'time':
   select_time_slice=slice(pd.to_datetime(check_sdate,format="%Y%m%d%H"),
                           pd.to_datetime(check_edate,format="%Y%m%d%H"))
   tmpds0=ds0.sel(time=select_time_slice).mean(dim=['lat','lon'])
else:
   tmpds0=ds0.mean(dim=['lat','lon'])

pltdf0=tmpds0.to_dataframe()


fig,ax=plt.subplots()
set_size(axe_w,axe_h,b=0.15)
#ax.set_prop_cycle(linestyle=['-','--','--'])
pltdf0[['btd_mean','btd_max','btd_min']].plot(ax=ax,marker='o',alpha=0.8,ms=2,zorder=4)
ax.legend(['Mean','Max','Min'])
ax.set_xlabel('%s [%s]' %(pltdf0.index.name.capitalize(),'$\mathrm{cm^{-1}}$'))
ax.set_ylabel('BT Differences [%s]' %(units))
xmin,xmax=ax.get_xbound()
ax.hlines(0.,xmin,xmax,colors='k',lw=0.8,ls='--',zorder=3)
ax.set_xlim(xmin,xmax)
#ax.hlines(0.,0,1,transform=ax.get_yaxis_transform(),colors='k',lw=0.8,ls='--',zorder=3)

outname='%s/%s_%s_%s.%s_%s.%s' %(imgsavpath,expn,sensor,pltvar,check_sdate,check_edate,ffmt)
if (fsave): print(outname,flush=1) ; fig.savefig(outname,dpi=quality); plt.close()

