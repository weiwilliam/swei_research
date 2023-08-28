#!/usr/bin/env python3
import os,sys
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mpcrs
import cartopy.crs as ccrs
import setuparea as setarea
from plot_utils import setupax_2dmap, plt_x2y, set_size

# Projection setting
proj=ccrs.PlateCarree(globe=None)

idxrg = range(9500,9600)

infile = '/data/users/swei/Dataset/jedi-obs/tmp/TROPOMI_S5P_20210815T18_CO_total.nc'
outfile = '/data/users/swei/FTPdir/test.png'

area='Glb'
minlon, maxlon, minlat, maxlat, crosszero, cyclic=setarea.setarea(area)
print(minlat,maxlat,minlon,maxlon,crosszero,cyclic)
cornll=[minlat,maxlat,minlon,maxlon]

ds = nc.Dataset(infile)

lat = ds.groups['MetaData'].variables['latitude'][:].ravel()
lon = ds.groups['MetaData'].variables['longitude'][:].ravel()

fig,ax,gl=setupax_2dmap(cornll,area,proj,12)
set_size(6,3)
sc=ax.scatter(lon[idxrg],lat[idxrg],c='b',s=2)

fig.savefig(outfile,dpi=300)
plt.close()

