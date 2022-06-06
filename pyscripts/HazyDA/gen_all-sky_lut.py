# -*- coding: utf-8 -*-
"""
Created on Tue Jun 29 16:06:23 2021

@author: ck102

Aerosol detection based on CADS 3.1 from NWP SAF

"""
import sys, os, platform
machine='S4'
os_name=platform.system()
if (machine=='MBP'):
    rootpath='/Users/weiwilliam'
    rootarch='/Volumes/WD2TB/ResearchData'
elif (machine=='Desktop'):
    rootpath='F:\GoogleDrive_NCU\Albany'
    rootarch='F:\ResearchData'
    rootgit='F:\GitHub\swei_research'
elif (machine=='S4'):
    rootpath='/data/users/swei'
    rootarch='/data/users/swei/ResearchData'
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
import numpy as np
import xarray as xa
import pandas as pd

# Plotting setup
sensor='iasi_metop-a'
lutfmt='csv'
ver='v5'
nchs=616

# Data path setup
lutpath=rootpath+'/AlbanyWork/Prospectus/Experiments/AeroObsStats/SD_LUT'

satstats_file=lutpath+'/'+sensor+'_'+str(nchs)+'_stats.'+ver+'.'+lutfmt
if (lutfmt=='xlsx'):
   lutdf=pd.read_excel(satstats_file)
elif (lutfmt=='csv'):
   lutdf=pd.read_csv(satstats_file)

filter = (lutdf.Aer_sen==1.)&(lutdf.iuse==1.)

nrows=lutdf.shape[0]
#lutdf['ich']=np.arange(1,nrows+1)

tmpdf=lutdf[['ich','Aeff_1','Aeff_2','SD_max']]
tmpdf=tmpdf.loc[filter,:]
#tmpdf=tmpdf.set_index('nuchan')
outputfile=lutpath+'/Ae_dep_SD.'+ver+'.txt'
np.savetxt(outputfile,tmpdf.values,fmt='%5i  %6.3f  %6.3f  %6.3f')
        
