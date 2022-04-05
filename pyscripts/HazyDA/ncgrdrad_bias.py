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
import setuparea as setarea
import numpy as np
from utils import ndate
from datetime import datetime
from datetime import timedelta
import matplotlib.dates as mdates
import xarray as xa
import pandas as pd
import scipy.stats as ss

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

biasterm=7
biastermname=['BC_Total',
              'BC_angord',
              'BC_Cloud_Liquid_Water',
              'BC_Constant',
              'BC_Cosine_Latitude_times_Node',
              'BC_Emissivity',
              'BC_Fixed_Scan_Position',
              'BC_Lapse_Rate',
              'BC_Lapse_Rate_Squared',
              'BC_Scan_Angle',
              'BC_Sine_Latitude']
#degres=2.5
degres=1

expset=1
if (expset==1):
   explist=np.array(['hazyda_ctrl','hazyda_aero'])
   leglist=['CTL','AER']
elif (expset==2):
   explist=np.array(['prctrl','praero'])
   leglist=['CTL','CAER']

sensor='iasi_metop-a'

sdate=2020061000
edate=2020071018
hint=6
chkwvn=962.5

area='Glb'
minlon, maxlon, minlat, maxlat, crosszero, cyclic=setarea.setarea(area)
print(minlat,maxlat,minlon,maxlon,crosszero,cyclic)

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
    
inpath=rootpath+'/archive/HazyDA/gridded_diag'

grdfile0='%s/%s_%s_%s_%s_bcterm%s_%.1fx%.1f.%s_%s.nc' %(inpath,leglist[0],sensor,loop,qcflg,biasterm,degres,degres,sdate,edate)
grdfile1='%s/%s_%s_%s_%s_bcterm%s_%.1fx%.1f.%s_%s.nc' %(inpath,leglist[1],sensor,loop,qcflg,biasterm,degres,degres,sdate,edate)

ds0=xa.open_dataset(grdfile0)
ds1=xa.open_dataset(grdfile1)


