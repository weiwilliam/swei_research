# -*- coding: utf-8 -*-
"""
Created on Tue Jun 29 16:06:23 2021

@author: ck102

Aerosol detection based on CADS 3.1 from NWP SAF

"""
import sys, os, platform
import numpy as np
import xarray as xa
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mpcrs
import cartopy.crs as ccrs
os_name=platform.system()
if (os_name=='Darwin'):
    rootpath='/Users/weiwilliam'
    rootarch='/Volumes/WD2TB/ResearchData'
elif (os_name=='Windows'):
    rootpath='F:\GoogleDrive_NCU\Albany'
    rootarch='F:\ResearchData'
    rootgit='F:\GitHub\swei_research'
elif (os_name=='Linux'):
    rootpath='/scratch2/BMC/gsd-fv3-dev/Shih-wei.Wei'
    rootarch='/scratch2/BMC/gsd-fv3-dev/Shih-wei.Wei/ResearchData'
    rootgit='/home/Shih-wei.Wei/research'
sys.path.append(rootgit+'/pyscripts/functions')
import setuparea as setarea
from plot_utils import setupax_2dmap, plt_x2y, set_size
from utils import ndate,setup_cmap
from datetime import datetime, timedelta
import scipy.stats

tlsize=12 ; lbsize=10
mpl.rc('axes', titlesize=tlsize,labelsize=lbsize)
mpl.rc('xtick',labelsize=lbsize)
mpl.rc('ytick',labelsize=lbsize)
mpl.rc('legend',fontsize='large')
fsave=1 ; ffmt='png' ; ptsize=4
axe_w=3 ; axe_h=3 ; quality=300

# Projection setting
proj=ccrs.PlateCarree(globe=None)

# Plotting setup
exp='AerObserver'
sensor='iasi_metop-a'
spectral_range=slice(700,1300)
loop='ges' #ges,anl
#if loop=='anl':
#    tlstr='OMA'
#elif loop=='ges':
#    tlstr='OMF'
plthist=0 # plot pdf of normalized OMB
plthist_mean_sd=1 # plot 2d histogram
lutfmt='csv'

# Data path setup
archpath=rootarch+'/Prospectus/AeroObsStats/nc_diag'
lutpath=rootpath+'/AlbanyWork/Prospectus/Experiments/AeroObsStats/AerObserver_newQC/SD_LUT/All'
outpath=rootpath+'/AlbanyWork/Prospectus/Experiments/AeroObsStats/images/'
archdir=archpath+'/'+exp

satstats_file=lutpath+'/'+sensor+'_'+str(nchs)+'_stats.'+lutfmt
if (lutfmt=='xlsx'):
   lutdf=pd.read_excel(satstats_file)
elif (lutfmt=='csv'):
   lutdf=pd.read_csv(satstats_file)

filter = (df.SD_o<df.SD_max)&(df.iuse==1.)

nrows=lutdf.shape[0]
df['ich']=np.arange(1,nrows+1)

tmpdf=lutdf[['ich','Aeff_1','Aeff_2','SD_max']]
tmpdf=tmpdf.loc[filter,:]
#tmpdf=tmpdf.set_index('nuchan')
np.savetxt('./test.txt',tmpdf.values,fmt='%5i  %6.3f  %6.3f  %6.3f')
        
