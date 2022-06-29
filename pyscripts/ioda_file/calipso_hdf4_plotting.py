# few different points in python to fortran
# 1. array indices start from 0
# 2. get specific value in an array by [] instead of ()
# 3. no ending statement for for-loop and if-statement
# 4. use [] for list variable
#
# import the necessary library,
# import [library].[module] as [shortname]
# from [library] import [module]
# then, can use the module with shortname
#
# hdf5 library, http://docs.h5py.org/en/stable/
import h5py as h5

# plotting tool, https://matplotlib.org/index.html
# available colormaps in matplotlib https://matplotlib.org/examples/color/colormaps_reference.html
import matplotlib.pyplot as plt 

# numpy package, https://numpy.org/doc/stable/ 
import numpy as np
 
# package to plot overlay on map, https://scitools.org.uk/cartopy/docs/latest/
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

#inpath='/data/users/clu/OBS_for_IODA/CALIPSO'
# Setup datapath and filename
inpath='/data/users/swei/Dataset/CALIPSO'
filename='CAL_LID_L2_05kmAPro-Standard-V4-20.2019-07-21T23-36-45ZD.h5'

pltvar='Backscatter_Coefficient_1064'
cmapname='gist_ncar'
nxticks=6 ; nyticks=10

# generate equal space array with 100 numbers between vmin and vmax
vmin=1e-4 ; vmax=1e-1
cnlvs=np.linspace(vmin,vmax,100) 

# load hdf5 data to read ('r')
# ipython: check available variables with h5ds.keys() 
h5ds=h5.File(inpath+'/'+filename,'r')

# pull metadata in hdf5 file
meta=h5ds['metadata']

# get the data altitude from metadata
height=meta['Lidar_Data_Altitudes'][0,:]

# get latitude and longitude
#   0: average of 15 shots
#   1: mid-point (8th) of 15 shots
#   2: point for final pulse
lat=h5ds['Latitude'][:,1]
lon=h5ds['Longitude'][:,1]

# plotting track
proj=ccrs.PlateCarree()    # set the projection
# Can use 0=False, 1=True for logical statement
fig=plt.figure(figsize=(8,5),constrained_layout=1)
ax=plt.subplot(projection=proj)
ax.coastlines(resolution='110m')
ax.set_extent((-180,180,-90,90),crs=proj)
ax.set_yticks([-90.,-60.,-30.,0.,30.,60.,90.],crs=proj)
ax.set_xticks([-180.,-120.,-60.,0.,60.,120.,180.],crs=proj)
ax.xaxis.set_major_formatter(LongitudeFormatter())
ax.yaxis.set_major_formatter(LatitudeFormatter())
ax.scatter(lon,lat,s=1,c='r')
fig.savefig('calipso_track.png',dpi=300)
plt.close()

#
# save the data in rawdata and swap dimension from (t,z) to (z,t)
# because contourf dimension has to be y-axis then x-axis
# np.where(cond,x,y), if cond is true, set value to x. cond is false set value to y
# 
rawdata=h5ds[pltvar][()].swapaxes(0,1)
mask=(rawdata!=-9999.0)
mskdata=np.where(mask,rawdata,np.nan)

# Set labels and ticks position of x- and y-yaxis 
xaxis=np.arange(mskdata.shape[1])
xtickspos=xaxis[::int(xaxis.size/(nxticks-1))]  
xlabels=[]
for x in xtickspos:
    if (lat[x]>0):
       ns='ºN'
    elif (lat[x]<0):
       ns='ºS'
    else:
       ns='º'
    if (lon[x]>0):
       ew='ºE'
    elif (lat[x]<0):
       ew='ºW'
    else:
       ew='º'
    xlabels.append('Lat %.1f%s\nLon %.1f%s' %(abs(lat[x]),ns,abs(lon[x]),ew))

yaxis=np.arange(height.size)
ytickspos=yaxis[::int(height.size/(nyticks-1))]
ylabels=[]
for y in ytickspos:
    ylabels.append('%.1f' %(height[y]))

r_asp=(5/height.size)/(9/xaxis.size)

# Plot with imshow
fig=plt.figure(figsize=(9,5),constrained_layout=1)
ax=plt.subplot()
# imshow
# https://matplotlib.org/api/_as_gen/matplotlib.axes.Axes.imshow.html?highlight=imshow#matplotlib.axes.Axes.imshow
# default aspect ratio is 1:1, need to adjust for better visualization
im=ax.imshow(mskdata,vmin=vmin,vmax=vmax,cmap=cmapname,aspect=r_asp)
ax.set_xticks(xtickspos)
ax.set_xticklabels(xlabels)
ax.set_yticks(ytickspos)
ax.set_yticklabels(ylabels)
ax.set_ylabel('Height[km]')
ax.set_title(pltvar,loc='left')
fig.colorbar(im,ax=ax,orientation='vertical',fraction=0.03)
fig.savefig('%s_imshow.png'%(pltvar),dpi=300)
plt.close()

# Plot with contourf
fig=plt.figure(figsize=(9,5),constrained_layout=1)
ax=plt.subplot()
# contourf(X,Y,Z[Y,X],**kwargs)
# https://matplotlib.org/api/_as_gen/matplotlib.axes.Axes.contourf.html?highlight=contourf#matplotlib.axes.Axes.contourf
cn=ax.contourf(xaxis,height,mskdata,levels=cnlvs,cmap=cmapname)
ax.set_xticks(xtickspos)
ax.set_xticklabels(xlabels)
ax.set_ylabel('Height[km]')
ax.set_title(pltvar,loc='left')
fig.colorbar(cn,ax=ax,orientation='vertical',fraction=0.03)
fig.savefig('%s_contourf.png'%(pltvar),dpi=300)
plt.close()
