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
elif (machine=='S4'):
    rootpath='/data/users/swei/AlbanyWork/Prospectus/Experiments/HazyDA/Images'
    rootarch='/data/users/swei'
elif (machine=='Hera'):
    rootpath='/scratch2/BMC/gsd-fv3-dev/Shih-wei.Wei'
    rootarch='/scratch2/BMC/gsd-fv3-dev/Shih-wei.Wei/ResearchData'
elif (machine=='Cheyenne'):
    rootpath='/glade/work/swei/output/images'
    rootarch='/scratch2/BMC/gsd-fv3-dev/Shih-wei.Wei/ResearchData'
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mpcrs
import xarray as xa
import pandas as pd

import setuparea as setarea
from plot_utils import plt_x2y, set_size
from utils import setup_cmap,find_cnlvs

tlsize=10 ; txsize=10
mpl.rc('axes', titlesize=tlsize,labelsize=txsize)
mpl.rc('xtick',labelsize=txsize)
mpl.rc('ytick',labelsize=txsize)
mpl.rc('legend',fontsize='small')
fsave=1 ; ffmt='png' ; ptsize=4
axe_w=6 ; axe_h=3 ; quality=300
tkfreq=2
minussign=u'\u2212'

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

pltvar='omb_nbc'
if pltvar=='omb_bc':
   varlb='OMB w/ BC'
elif pltvar=='omb_nbc':
   varlb='OMB w/o BC'
sdate=2020060106
edate=2020071018
hint=6

sensor='iasi_metop-a'
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
outputpath=rootpath+'/DiagFiles/gridded_rad'
imgsavpath=outputpath+'/spectrum/omb/'+area
if ( not os.path.exists(imgsavpath) ):
   os.makedirs(imgsavpath)

grdfile0='%s/%s_%s_%s_%s_omb_%.1fx%.1f.mean.%s_%s.nc' %(inpath,expnlist[0],sensor,loop,qcflg,degres,degres,sdate,edate)
grdfile1='%s/%s_%s_%s_%s_omb_%.1fx%.1f.mean.%s_%s.nc' %(inpath,expnlist[1],sensor,loop,qcflg,degres,degres,sdate,edate)

ds0=xa.open_dataset(grdfile0)
ds1=xa.open_dataset(grdfile1)

sum_dims=['lat','lon']

cnt_name='%s_%s'%(pltvar,'count')
mean_name='%s_%s'%(pltvar,'mean')
var_name='%s_%s'%(pltvar,'var')
cnt0=ds0[cnt_name]; cnt1=ds1[cnt_name]

mean0=(ds0[mean_name]*cnt0).sum(dim=sum_dims)/(cnt0.sum(dim=sum_dims))
mean1=(ds1[mean_name]*cnt1).sum(dim=sum_dims)/(cnt1.sum(dim=sum_dims))

var0=(ds0[var_name]*cnt0).sum(dim=sum_dims)/(cnt0.sum(dim=sum_dims))
var1=(ds1[var_name]*cnt1).sum(dim=sum_dims)/(cnt1.sum(dim=sum_dims))

rmse0=np.sqrt(var0+mean0*mean0)
rmse1=np.sqrt(var1+mean1*mean1)

wvn=ds0.wavenumber.data
wvl=1e4/wvn
wvnlb='Wavenumber [$cm^{-1}$]'
wvllb='Wavelength [µm]'
prop_dict={'color'     :['b','r'],
           'line_style':['-','-'],
           'line_width':[1.5,1.5],
           'marker'    :['o','o'],
           'mark_size' :[3.,3.],
           'legend'    :expnlist,
           }
tistr=''

for stat_type in ['Mean','RMS']:
    if stat_type=='Mean':
       tmpda=xa.concat((mean0,mean1),dim='exps')
    elif stat_type=='RMS':
       tmpda=xa.concat((rmse0,rmse1),dim='exps')
    pltda=tmpda.data.swapaxes(0,1)
    yaxlb='%s %s [%s]' %(stat_type,varlb,units)
    
    fig,ax=plt.subplots()
    set_size(axe_w,axe_h,ax=ax,b=0.25)
    plt_x2y(pltda,yaxlb,wvn,wvnlb,wvl,wvllb,prop_dict,tistr,0,[],ax=ax)
    
    if (fsave):
       outname='%s/%s_%s-%s_%s_%s.%s_%s.%s' %(imgsavpath,stat_type,expnlist[1],expnlist[0],sensor,pltvar,sdate,edate,ffmt)
       print(outname,flush=1)
       fig.savefig(outname,dpi=quality)
       plt.close()


