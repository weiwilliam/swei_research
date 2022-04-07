# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 21:03:35 2022

@author: ck102
"""
import numpy as np
import xarray as xa
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mpcrs
from matplotlib.gridspec import GridSpec
import cartopy.crs as ccrs
import cartopy.feature as cft
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from scipy.interpolate import interp1d
from netCDF4 import Dataset

# label size of plot
lbsize=11; quality=500
hcbori='horizontal'; hcb_frac=0.03;  hcb_pad=0.08
vcbori='vertical'  ; vcb_frac=0.025; vcb_pad=0.06


# Full Disk, CONUS, Regional
filepath="C:/Users/ck102/Documents/Globus/GOES-16"
filename=filepath+'/OR_ABI-L2-ADPC-M6_G16_s20220892351172_e20220892353545_c20220892355250.nc'

readin=Dataset(filename)
### Convert radian to lat/lon
# Degree conversion code is from the webpage below
# https://www.star.nesdis.noaa.gov/smcd/spb/aq/AMS_Short_Course/abi_visualize.php
# Ignore numpy error for sqrt of negative number ('x'); 
# occurs for GOES-16 ABI CONUS sector data
np.seterr(invalid='ignore')

# Read in GOES Imager Projection data
lat_rad_1d = readin['x'][:]
lon_rad_1d = readin['y'][:]
projection_info = readin['goes_imager_projection']
lon_origin = projection_info.longitude_of_projection_origin
H = projection_info.perspective_point_height+projection_info.semi_major_axis
r_eq = projection_info.semi_major_axis
r_pol = projection_info.semi_minor_axis

# Create meshgrid filled with radian angles
lat_rad,lon_rad = np.meshgrid(lat_rad_1d,lon_rad_1d)

# lat/lon calculus routine from satellite radian angle vectors
lambda_0 = (lon_origin*np.pi)/180.0

a_var = np.power(np.sin(lat_rad),2.0) + (np.power(np.cos(lat_rad),2.0)*(np.power(np.cos(lon_rad),2.0)+(((r_eq*r_eq)/(r_pol*r_pol))*np.power(np.sin(lon_rad),2.0))))
b_var = -2.0*H*np.cos(lat_rad)*np.cos(lon_rad)
c_var = (H**2.0)-(r_eq**2.0)

r_s = (-1.0*b_var - np.sqrt((b_var**2)-(4.0*a_var*c_var)))/(2.0*a_var)

s_x = r_s*np.cos(lat_rad)*np.cos(lon_rad)
s_y = - r_s*np.sin(lat_rad)
s_z = r_s*np.cos(lat_rad)*np.sin(lon_rad)

lat = (180.0/np.pi)*(np.arctan(((r_eq*r_eq)/(r_pol*r_pol))*((s_z/np.sqrt(((H-s_x)*(H-s_x))+(s_y*s_y))))))
lon = (lambda_0 - np.arctan(s_y/(H-s_x)))*(180.0/np.pi)

lon_bound = [-152.10928, -52.94688]
lat_bound = [14.57134, 56.76145]
x_bound = ds['x_image_bounds'][:]
y_bound = ds['y_image_bounds'][:]
x_interp = interp1d(x_bound, lon_bound)
y_interp = interp1d(y_bound, lat_bound)
alon = x_interp(np.array(ds['x'][:]))
alat = y_interp(np.array(ds['y'][:]))
aalon, aalat = np.meshgrid(alon,alat)

# Read in dust and smoke mask
smoke_msk=ds.Smoke
dust_msk=ds.Dust
qualityflag=ds.DQF
#
smoke_arr=np.zeros_like(lat)
low_c_smoke=(smoke_msk==1)&(qualityflag==0)
med_c_smoke=(smoke_msk==1)&(qualityflag==4)
hi_c_smoke= (smoke_msk==1)&(qualityflag==12)
smoke_arr=np.where(low_c_smoke,1,smoke_arr)
smoke_arr=np.where(med_c_smoke,2,smoke_arr)
smoke_arr=np.where(hi_c_smoke, 3,smoke_arr)
pltsmk=(low_c_smoke)|(med_c_smoke)|(hi_c_smoke)

dust_arr=np.zeros_like(lat)
low_c_dust=(dust_msk==1)&(qualityflag==0)
med_c_dust=(dust_msk==1)&(qualityflag==16)
hi_c_dust= (dust_msk==1)&(qualityflag==48)
dust_arr=np.where(low_c_dust,1,dust_arr)
dust_arr=np.where(med_c_dust,2,dust_arr)
dust_arr=np.where(hi_c_dust, 3,dust_arr)
pltdst=(low_c_dust)|(med_c_dust)|(hi_c_dust)
#
aer_conf_lvs=[0.5,1.5,2.5,3.5]
aer_norm = mpcrs.BoundaryNorm(aer_conf_lvs,4)
smkcmap=mpcrs.ListedColormap(['pink','deeppink','red'])
dstcmap=mpcrs.ListedColormap(['khaki','peru','saddlebrown'])
cbtck_lbs=['Low','Med','High']
smkcbn='Smoke Confidence'
dstcbn='Dust Confidence'

# Title
timestr=pd.to_datetime(ds.time_bounds.mean().data).strftime('%Y-%m-%d %H:%M:%SZ')
titlestr='%s ADP Smoke/Dust Confidence Flag' %(timestr)

# Setup projections
proj_pc=ccrs.PlateCarree()
proj=ccrs.Geostationary(satellite_height=projection_info.perspective_point_height,
                        central_longitude=projection_info.longitude_of_projection_origin,
                        sweep_axis=projection_info.sweep_angle_axis)

# Generate full-disk image
fig=plt.figure(figsize=(6,8))
gs=GridSpec(2,2,width_ratios=(1,1),height_ratios=(25,1),hspace=0.)
ax=fig.add_subplot(gs[:-1,:],projection=proj)
cax1=fig.add_subplot(gs[1,0])
cax2=fig.add_subplot(gs[1,1])
ax.coastlines(resolution='110m',color='cyan')
ax.set_global()
ax.set_facecolor('k')
ax.add_feature(cft.STATES,edgecolor='cyan')
ax.set_title(titlestr)
gl=ax.gridlines(draw_labels=True,dms=True,x_inline=True, y_inline=False)
gl.right_labels=False
gl.top_labels=False
gl.xformatter=LongitudeFormatter(degree_symbol=u'\u00B0 ')
gl.yformatter=LatitudeFormatter(degree_symbol=u'\u00B0 ')
gl.xlabel_style={'size':lbsize,'color':'white'}
gl.ylabel_style={'size':lbsize}

sc1=ax.scatter(lon[pltsmk==1],lat[pltsmk==1],s=0.2,
            c=smoke_arr[pltsmk==1],cmap=smkcmap,
            norm=aer_norm,zorder=3,transform=proj_pc)
cb1=plt.colorbar(sc1,cax=cax1,ticks=[1.,2.,3.],label=smkcbn,
                 orientation=hcbori)
cb1.ax.set_xticklabels(cbtck_lbs)

sc2=ax.scatter(lon[pltdst==1],lat[pltdst==1],s=0.2,
            c=smoke_arr[pltdst==1],cmap=dstcmap,
            norm=aer_norm,zorder=3,transform=proj_pc)
cb2=plt.colorbar(sc2,cax=cax2,ticks=[1.,2.,3.],label=dstcbn,
                 orientation=hcbori)
cb2.ax.set_xticklabels(cbtck_lbs)

# Define the output name and save
imgname=filepath+'/FD_test.png'
fig.savefig(imgname,dpi=quality)
plt.close()

# Generate the plot for regional
minlon=np.nanmin(lon)-0.2; maxlon=np.nanmax(lon)+0.2
minlat=np.nanmin(lat)-0.2; maxlat=np.nanmax(lat)+0.2

fig=plt.figure(figsize=(7.5,4.5))
gs=GridSpec(2,2,width_ratios=(1,1),height_ratios=(20,1),hspace=0.)
ax=fig.add_subplot(gs[:-1,:],projection=proj_pc)
cax1=fig.add_subplot(gs[1,0])
cax2=fig.add_subplot(gs[1,1])
ax.coastlines(resolution='110m',color='cyan')
ax.set_extent((minlon,maxlon,minlat,maxlat),crs=proj_pc)
ax.set_facecolor('k')
ax.add_feature(cft.STATES,edgecolor='cyan')
ax.set_title(titlestr)
gl=ax.gridlines(draw_labels=True,dms=True,x_inline=False, y_inline=False)
gl.right_labels=False
gl.top_labels=False
gl.xformatter=LongitudeFormatter(degree_symbol=u'\u00B0 ')
gl.yformatter=LatitudeFormatter(degree_symbol=u'\u00B0 ')
gl.xlabel_style={'size':lbsize}
gl.ylabel_style={'size':lbsize}

sc1=ax.scatter(aalon[pltsmk==1],aalat[pltsmk==1],s=0.2,
            c=smoke_arr[pltsmk==1],cmap=smkcmap,
            norm=aer_norm,zorder=3)
cb1=plt.colorbar(sc1,cax=cax1,ticks=[1.,2.,3.],label=smkcbn,
                 orientation=hcbori)
cb1.ax.set_xticklabels(cbtck_lbs)

sc2=ax.scatter(aalon[pltdst==1],aalat[pltdst==1],s=0.2,
            c=smoke_arr[pltdst==1],cmap=dstcmap,
            norm=aer_norm,zorder=3)
cb2=plt.colorbar(sc2,cax=cax2,ticks=[1.,2.,3.],label=dstcbn,
                 orientation=hcbori)
cb2.ax.set_xticklabels(cbtck_lbs)

# Define the output name and save
imgname=filepath+'/Reg_test1.png'
fig.savefig(imgname,dpi=quality)
plt.close()




