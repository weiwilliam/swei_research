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
from matplotlib.dates import (DAILY, DateFormatter,
                              rrulewrapper, RRuleLocator)

# Plotting setup
mpl.rc('axes',titlesize=14,labelsize=12)
mpl.rc('xtick',labelsize=12)
mpl.rc('ytick',labelsize=12)
mpl.rc('legend',fontsize='medium')
axe_w=6; axe_h=3
quality=300
tkfreq=1
clrlst=['black','tab:red','tab:blue']
lstylst=['-','-','-']
mrkrlst=[' ',' ',' ']

inputpath='/data/users/swei/archive'
outputpath=rootpath+'/Ens_size/Timeseries/RMSEvsSPRD'
if ( not os.path.exists(outputpath) ):
    os.makedirs(outputpath)

sdate=2020060300
edate=2020061218
hint=6
pltvar='total'
area='Glb'
explist=['ctrl','10ens','30ens']
levlist=[850.,500.,250.]
varlist=['tmp','spfh','ugrd','vgrd']

#
cnlvs=np.array((0., 0.05, 0.1, 0.15, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 2.5))
clridx=np.array((0,2,4,7,8,9,10,12,14,15,16,17,18))
clrmap=setup_cmap('precip3_16lev',clridx)
aer_norm = mpcrs.BoundaryNorm(cnlvs,len(clridx),extend='both')
#cblb='RMSE/Spread'

# Constant configuration
proj=ccrs.PlateCarree()
grav=9.80665e+0

# Setup plotting area
minlon, maxlon, minlat, maxlat, crosszero, cyclic=setarea.setarea(area)
print(minlat,maxlat,minlon,maxlon,crosszero,cyclic)
if (area=='Glb'):
   minlon=-180. ; maxlon=180.
cornerll=[minlat,maxlat,minlon,maxlon]
# 
syy=int(str(sdate)[:4]); smm=int(str(sdate)[4:6])
sdd=int(str(sdate)[6:8]); shh=int(str(sdate)[8:10])
eyy=int(str(edate)[:4]); emm=int(str(edate)[4:6])
edd=int(str(edate)[6:8]); ehh=int(str(edate)[8:10])

date1 = datetime(syy,smm,sdd,shh)
date2 = datetime(eyy,emm,edd,ehh)
delta = timedelta(hours=hint)
dates = pd.date_range(start=date1, end=date2, freq=delta)

tnum=0
dlist=[]
cdate=sdate
while (cdate<=edate):
    dlist.append(str(cdate))
    tnum=tnum+1
    cdate=ndate(hint,cdate)

rule = rrulewrapper(DAILY, byhour=6, interval=2)
loc = RRuleLocator(rule)
formatter = DateFormatter('%Y%h %n %d %Hz')

for var in varlist:
    for lev in levlist:
        for date in dlist:
            rmse_file=inputpath+'/Ens_size/rmse.'+date+'.nc4'
            sprd_file=inputpath+'/Ens_size/sprd.'+date+'.nc4'
            rmse=xa.open_dataset(rmse_file)
            sprd=xa.open_dataset(sprd_file)

            rmsetmp=rmse[var].sel(pfull=lev,method='nearest')
            sprdtmp=sprd[var].sel(pfull=lev,method='nearest')
            rmse_mean_tmp=rmsetmp.mean(dim=['grid_xt','grid_yt'])
            sprd_mean_tmp=sprdtmp.mean(dim=['grid_xt','grid_yt'])

            if (date==dlist[0]):
               rmse_mean=rmse_mean_tmp
               sprd_mean=sprd_mean_tmp
            else:
               rmse_mean=xa.concat((rmse_mean,rmse_mean_tmp),dim='time')
               sprd_mean=xa.concat((sprd_mean,sprd_mean_tmp),dim='time')
    
        outname='%s/TS_%s_%s_%s.png' %(outputpath,area,lev,var)
        pltrmse=rmse_mean.values.swapaxes(0,1)
        pltsprd=sprd_mean.values.swapaxes(0,1)
        pltdates=rmse_mean.time
        fig,ax=plt.subplots()
        set_size(axe_w,axe_h,b=0.13)
        ax.set_prop_cycle(color=clrlst) #,linestyle=lstylst,marker=mrkrlst)
        ax.xaxis.set_major_locator(loc)
        ax.xaxis.set_major_formatter(formatter)
        ax.plot_date(pltdates,pltrmse,'-')
        ax.plot_date(pltdates,pltsprd,'--')
        lglst=[]
        for val in ['rmse','sprd']:
            for expn in rmse.exps.values:
                lglst.append('%s_%s'%(expn,val))
        ax.legend(lglst)
        if (var=='spfh'):
            ax.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
        ax.set_ylabel('%s hPa %s' %(lev,var))
        ax.grid()
        print(outname)
        fig.savefig(outname,dpi=quality)
        plt.close()

        outname='%s/TS_%s_%s_%s.ratio.png' %(outputpath,area,lev,var)
        pltdata=pltsprd/pltrmse
        fig,ax=plt.subplots()
        set_size(axe_w,axe_h,b=0.13)
        ax.set_prop_cycle(color=clrlst) #,linestyle=lstylst,marker=mrkrlst)
        ax.xaxis.set_major_locator(loc)
        ax.xaxis.set_major_formatter(formatter)
        ax.plot_date(pltdates,pltdata,'o-')
        ax.legend(rmse.exps.values)
        ax.set_ylabel('%s hPa %s' %(lev,var))
        ax.grid()
        print(outname)
        fig.savefig(outname,dpi=quality)
        plt.close()

