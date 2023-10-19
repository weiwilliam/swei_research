#!/usr/bin/env python3
import sys, os, platform
import numpy as np
import netCDF4 as nc
import pandas as pd
import matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mpcrs
from datetime import datetime, timedelta
from matplotlib.colors import LinearSegmentedColormap
from plot_utils import set_size, setupax_2dmap
import setuparea as setarea
import cartopy.crs as ccrs
#
axe_w=6; axe_h=3;
quality=300; ffmt='png'
# Projection setting
proj=ccrs.PlateCarree(globe=None)

cdate=2016070100
cymstr=str(cdate)[:6]

explist=['noda','anal']
obtype='pm25'

outpath='/data/users/swei/MAPP/pm25_scimg'
archdir='/data/users/swei/MAPP/pm25_hofx'
savedir=outpath+'/pm25_scatter'
suffix='nc4'

def hofx_dict(ncd,varname):
    outdict = {'nlocs':  ncd.dimensions['Location'].size,
               'lats':   ncd.groups['MetaData'].variables['latitude' ][:],
               'lons':   ncd.groups['MetaData'].variables['longitude'][:],
               'obs':    ncd.groups['ObsValue'].variables[varname][:],
               'hofx':   ncd.groups['hofx'].variables[varname][:],
               'preqc':  ncd.groups['PreQC'].variables[varname][:],
               'effqc':  ncd.groups['EffectiveQC'].variables[varname][:],
               }
    return outdict

def process_pm25(in_dict):
    effqc_mask = in_dict['effqc'].data[:] == 0
    pm25_obs = in_dict['obs'].data[effqc_mask]
    pm25_hofx = in_dict['hofx'].data[effqc_mask]
    pm25_lats = in_dict['lats'].data[effqc_mask]
    pm25_lons = in_dict['lons'].data[effqc_mask]

    return pm25_obs, pm25_hofx, effqc_mask, pm25_lats, pm25_lons

area='Glb'
minlon, maxlon, minlat, maxlat, crosszero, cyclic=setarea.setarea(area)
print(minlat,maxlat,minlon,maxlon,crosszero,cyclic)
if (area=='Glb'):
   minlon=-180. ; maxlon=180.
cornll=[minlat,maxlat,minlon,maxlon]

file1 = ('%s/%s_%s/%s/openaq_pm25_hofx.%s.%s' % (archdir,explist[0],obtype,cymstr,cdate,suffix))
print('Processing: %s'%file1,flush=1)
ncd1 = nc.Dataset(file1)
dict1 = hofx_dict(ncd1,obtype)

pm25_obs1, pm25_hofx1, mask1, lats, lons = process_pm25(dict1)
print(lats.size)

fig,ax,gl=setupax_2dmap(cornll,area,proj,12)
set_size(axe_w,axe_h)
#ax.set_title(tistr,loc='left')
sc=ax.scatter(lons,lats,c='b',s=1.6)

fname=('./%s_openaq.%s.%s' %(area,str(cdate),ffmt))
fig.savefig(fname,dpi=quality)
plt.close()
