# -*- coding: utf-8 -*-
"""
Created on Tue Jun 29 16:06:23 2021

@author: ck102

Aerosol detection based on CADS 3.1 from NWP SAF

"""
import sys, os, platform
os_name=platform.system()
if (os_name=='Darwin'):
    rootpath='/Users/weiwilliam'
    rootarch='/Volumes/WD2TB/ResearchData'
elif (os_name=='Windows'):
    rootpath='F:\GoogleDrive_NCU\Albany'
    rootarch='F:\ResearchData'
    rootgit='F:\GitHub\swei_research'
elif (os_name=='Linux'):
    if (os.path.exists('/scratch1')):
        rootpath='/scratch2/BMC/gsd-fv3-dev/Shih-wei.Wei'
        rootarch='/scratch2/BMC/gsd-fv3-dev/Shih-wei.Wei/ResearchData'
        rootgit='/home/Shih-wei.Wei/research'
    elif (os.path.exists('/glade')):
        rootpath='/glade/work/swei/output/images'
        rootarch='/scratch2/BMC/gsd-fv3-dev/Shih-wei.Wei/ResearchData'
        rootgit='/glade/u/home/swei/research'
        machine='Cheyenne'
    elif (os.path.exists('/cardinal')):
        rootpath='/data/users/swei/Experiments/AeroObsStats'
        rootarch='/scratch/users/swei/ncdiag'
        rootgit='/home/swei/research'
        machine='S4'
sys.path.append(rootgit+'/pyscripts/functions')
import numpy as np
import xarray as xa
import pandas as pd

# Plotting setup
exp='AerObserver'
sensor='iasi_metop-a'
lutfmt='csv'
ver='v3'
nchs=616

# Data path setup
lutpath=rootpath+'/SD_LUT'

satstats_file=lutpath+'/'+sensor+'_'+str(nchs)+'_stats_new.'+ver+'.'+lutfmt
if (lutfmt=='xlsx'):
   lutdf=pd.read_excel(satstats_file)
elif (lutfmt=='csv'):
   lutdf=pd.read_csv(satstats_file)

filter = (lutdf.SD_o<lutdf.SD_max)&(lutdf.iuse==1.)

nrows=lutdf.shape[0]
lutdf['ich']=np.arange(1,nrows+1)

tmpdf=lutdf[['ich','Aeff_1','Aeff_2','SD_max']]
tmpdf=tmpdf.loc[filter,:]
#tmpdf=tmpdf.set_index('nuchan')
outputfile=lutpath+'/Ae_dep_SD.'+ver+'.txt'
np.savetxt(outputfile,tmpdf.values,fmt='%5i  %6.3f  %6.3f  %6.3f')
        
