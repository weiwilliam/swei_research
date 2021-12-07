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
    elif (os.path.exists('/cardinal')):
        rootpath='/data/users/swei/Images'
        rootarch='/scratch2/BMC/gsd-fv3-dev/Shih-wei.Wei/ResearchData'
        rootgit='/home/swei/research'
sys.path.append(rootgit+'/pyscripts/functions')
from utils import setup_cmap, ndate
from plot_utils import setupax_2dmap,set_size
import setuparea as setarea
import xarray as xa
import numpy as np
import pandas as pd
from datetime import datetime, timedelta
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mpcrs
import cartopy.crs as ccrs

# Plotting setup
mpl.rc('axes',titlesize=18,labelsize=16)
mpl.rc('xtick',labelsize=16)
mpl.rc('ytick',labelsize=16)
mpl.rc('legend',fontsize='xx-large')
axe_w=8; axe_h=3
quality=300
tkfreq=1
clrlst=['black','tab:red','tab:blue']
lstylst=['-','-','-']
mrkrlst=[' ',' ',' ']

inputpath='/data/users/swei/archive'
outputpath=rootpath+'/Ens_size/Timeseries'
if ( not os.path.exists(outputpath) ):
    os.makedirs(outputpath)

sdate=2020060106
edate=2020061218
hint=6
pltvar='total'
area='Glb'
explist=['ctrl','10ens','30ens']
levlist=[850.,500.]
varlist=['tmp','spfh','ugrd','vgrd']

#
cnlvs=np.array((0., 0.05, 0.1, 0.15, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 2.5))
clridx=np.array((0,2,4,7,8,9,10,12,14,15,16,17,18))
clrmap=setup_cmap('precip3_16lev',clridx)
aer_norm = mpcrs.BoundaryNorm(cnlvs,len(clridx),extend='both')
cblb='RMSE/Spread'

# Constant configuration
proj=ccrs.PlateCarree()
grav=9.80665e+0

# Setup plotting area
minlon, maxlon, minlat, maxlat, crosszero, cyclic=setarea.setarea(area)
print(minlat,maxlat,minlon,maxlon,crosszero,cyclic)
if (area=='Glb'):
   minlon=-180. ; maxlon=180.
cornerll=[minlat,maxlat,minlon,maxlon]

ens_file=inputpath+'/Ens_size/rmse_spread_ratio.nc4'
ratio=xa.open_dataset(ens_file)
dates=ratio.time

for var in varlist:
    for lev in levlist:
        tmp=ratio[var].sel(pfull=lev,method='nearest')
        tmp_std=tmp.std(dim=['grid_xt','grid_yt'])
        tmp_mean=tmp.mean(dim=['grid_xt','grid_yt'])

        outname='%s/TS_%s_%s_%s.png' %(outputpath,area,lev,var)
        pltdata=tmp_mean.values.swapaxes(0,1)
        fig,ax=plt.subplots()
        set_size(axe_w,axe_h)
        ax.set_prop_cycle(color=clrlst,linestyle=lstylst,marker=mrkrlst)
        ax.plot_date(dates,pltdata)
        ax.legend(ratio.exps.values)
        ax.set_title('%s %s' %(lev,var),loc='left')
        print(outname)
        fig.savefig(outname,dpi=quality)
        plt.close()
