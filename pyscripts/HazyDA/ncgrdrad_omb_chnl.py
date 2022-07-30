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

pltvar='omb_bc_mean'
varlb='OMB w/ BC'
sdate=2020060106
edate=2020071018
hint=6

sensor='iasi_metop-a'
chkwvn=962.5
degres=2.5
units='K'

explist=np.array(['hazyda_ctrl','hazyda_aero'])
expnlist=['CTL','AER']

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
outputpath=rootpath+'/DiagFiles/gridded'
imgsavpath=outputpath+'/2dmap/omb/'+area
if ( not os.path.exists(imgsavpath) ):
   os.makedirs(imgsavpath)

grdfile0='%s/%s_%s_%s_%s_omb_%.1fx%.1f.mean.%s_%s.nc' %(inpath,expnlist[0],sensor,loop,qcflg,degres,degres,sdate,edate)
grdfile1='%s/%s_%s_%s_%s_omb_%.1fx%.1f.mean.%s_%s.nc' %(inpath,expnlist[1],sensor,loop,qcflg,degres,degres,sdate,edate)

ds0=xa.open_dataset(grdfile0)
ds1=xa.open_dataset(grdfile1)

pltda0=ds0[pltvar].sel(wavenumber=chkwvn)
pltda1=ds1[pltvar].sel(wavenumber=chkwvn)

tmpda=xa.concat((pltda0,pltda1),dim='exps')

if (fsave):
   print(outname)
   fig.savefig(outname,dpi=quality)
   plt.close()


