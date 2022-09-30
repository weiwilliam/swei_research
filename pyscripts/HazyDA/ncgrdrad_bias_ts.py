#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 21:59:35 2019

@author: weiwilliam
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
sys.path.append(rootgit+'/pyscripts/functions')
import setuparea as setarea
import numpy as np
from datetime import datetime
from datetime import timedelta
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mpcrs
from matplotlib.dates import DateFormatter
import xarray as xa
import pandas as pd

import setuparea as setarea
from plot_utils import setupax_2dmap, plt_x2y, set_size
from utils import setup_cmap,find_cnlvs

tlsize=12 ; txsize=10
mpl.rc('axes', titlesize=tlsize,labelsize=txsize)
mpl.rc('xtick',labelsize=txsize)
mpl.rc('ytick',labelsize=txsize)
mpl.rc('legend',fontsize='small')
fsave=1 ; ffmt='png' ; ptsize=4
axe_w=6 ; axe_h=2.7 ; quality=300
tkfreq=2
minussign=u'\u2212'

pltcolor=['blue', 'red', 'black', 'grey']
pltstyle= ['-','-','-','-']
date_fmt= DateFormatter('%Y %h %n %d %Hz')


sensorlist=['airs_aqua','amsua_aqua','amsua_metop-a','amsua_n15','amsua_n18',
            'amsua_n19','atms_npp','avhrr_metop-a','avhrr_n18','cris_npp','gmi_gpm',
            'hirs4_metop-a','hirs4_metop-b','hirs4_n19','iasi_metop-a','iasi_metop-b',
            'mhs_metop-a','mhs_metop-b','mhs_n18','mhs_n19','saphir_meghat',
            'seviri_m08','seviri_m10','sndrd1_g15','sndrd2_g15','sndrd3_g15',
            'sndrd4_g15','ssmis_f17','ssmis_f18']
degres=2.5
#degres=1

explist=np.array(['hazyda_ctrl','hazyda_aero'])
expnlist=['CTL','AER']

sensor='iasi_metop-a'

sdate=2020060106
edate=2020071018
hint=6
chkwvn=962.5
nterm=6
pltvar='bcterm_mean'
units='K'

area='Glb'
minlon, maxlon, minlat, maxlat, crosszero, cyclic=setarea.setarea(area)
print(minlat,maxlat,minlon,maxlon,crosszero,cyclic)
if (area=='Glb'):
   minlon=-180. ; maxlon=180.
else:
   minlon=(minlon+180)%360-180
   maxlon=(maxlon+180)%360-180
cornerll=[minlat,maxlat,minlon,maxlon]

useqc=-1
if (useqc==-2):
    qcflg='noqc'
elif (useqc==-1):
    qcflg='qc'
else:
    qcflg='qc%s'%(useqc)

loop='ges' #ges,anl
if loop=='anl':
    tlstr='OMA'
elif loop=='ges':
    tlstr='OMF'

inpath=rootarch+'/archive/HazyDA/gridded_diag'
outputpath=rootpath+'/DiagFiles/gridded_rad'
imgsavpath=outputpath+'/1ch/bcterm'
if ( not os.path.exists(imgsavpath) ):
   os.makedirs(imgsavpath)

grdfile0='%s/%s_%s_%s_%s_bcterm_%.1fx%.1f.time.%s_%s.nc' %(inpath,expnlist[0],sensor,loop,qcflg,degres,degres,sdate,edate)
grdfile1='%s/%s_%s_%s_%s_bcterm_%.1fx%.1f.time.%s_%s.nc' %(inpath,expnlist[1],sensor,loop,qcflg,degres,degres,sdate,edate)

ds0=xa.open_dataset(grdfile0)
ds1=xa.open_dataset(grdfile1)
bctermname=ds0.nbcterm.data[nterm]

pltda0=ds0[pltvar].sel(nbcterm=bctermname,wavenumber=chkwvn)
pltda1=ds1[pltvar].sel(nbcterm=bctermname,wavenumber=chkwvn)

df0=pltda0.mean(dim=('lat','lon')).to_dataframe().rename(columns={pltvar:expnlist[0]})
df1=pltda1.mean(dim=('lat','lon')).to_dataframe().rename(columns={pltvar:expnlist[1]})

tmpdf=pd.concat((df0,df1),axis=1)[expnlist]

fig,ax=plt.subplots()
set_size(axe_w,axe_h,b=0.12)
ax.set_prop_cycle(color=pltcolor, linestyle=pltstyle)
tmpdf.plot(ax=ax,marker='o',ms=4)
ax.xaxis.set_major_formatter(date_fmt)
ax.grid(axis='x')
tistr='%s %.2f $\mathrm{cm^{-1}}$'%(sensor,chkwvn)
ax.set_title(tistr,loc='left')
ax.set_ylabel('%s [%s]'%(bctermname.replace('_',' '),units))
ax.set_xlabel('')

if (fsave):
   outname='%s/%s_%s_%s_%s_%.2f.png' %(imgsavpath,sensor,bctermname,expnlist[0],expnlist[1],chkwvn)
   print(outname,flush=1)
   fig.savefig(outname,dpi=quality)
   plt.close()
