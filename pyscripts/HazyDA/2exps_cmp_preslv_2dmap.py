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
    rootpath='/data/users/swei/AlbanyWork/Prospectus/Experiments/HazyDA/Images'
    rootarch='/data/users/swei'
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
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mpcrs
import matplotlib.dates as mdates
from matplotlib.dates import (DAILY, DateFormatter,
                              rrulewrapper, RRuleLocator)
import setuparea as setarea
from plot_utils import setupax_2dmap, plt_x2y, set_size
from utils import ndate,setup_cmap
from datetime import datetime, timedelta
import cartopy.crs as ccrs

tlsize=12 ; txsize=10
mpl.rc('axes', titlesize=tlsize,labelsize=txsize)
mpl.rc('xtick',labelsize=txsize)
mpl.rc('ytick',labelsize=txsize)
mpl.rc('legend',fontsize='small')
fsave=1 ; ffmt='png' ; ptsize=4
axe_w=6 ; axe_h=3 ; quality=300
tkfreq=2
minussign=u'\u2212'

# Projection setting
proj=ccrs.PlateCarree(globe=None)

outputpath=rootpath+'/grib'
expsarch=rootarch+'/archive'

explist=['hazyda_ctrl','hazyda_aero']
expnlist=['CTL','AER']
#enum=explist.shape[0]

sdate=2020061000
edate=2020071018
hint=6
pltvar='t' # z, r, q, t, u, v
units='K'  # m,'K','%','g/kg','K','m/s','mb'
grav=9.80665e+0

if (pltvar=='z'):
    expvar='gh'
else:
    expvar=pltvar

area='Glb'
minlon, maxlon, minlat, maxlat, crosszero, cyclic=setarea.setarea(area)
print(area,minlat,maxlat,minlon,maxlon,crosszero)
if (area=='Glb'):
   minlon=-180. ; maxlon=180.
else:
   minlon=(minlon+180)%360-180
   maxlon=(maxlon+180)%360-180
cornerll=[minlat,maxlat,minlon,maxlon]

imgsavpath=outputpath+'/2dmap/'+area
if ( not os.path.exists(imgsavpath) ):
   os.makedirs(imgsavpath)

cbori='horizontal' #vertical, horizontal
if (cbori=='vertical'):
   cb_frac=0.025
   cb_pad=0.06
elif (cbori=='horizontal'):
   cb_frac=0.04
   cb_pad=0.1

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

date_cnt=0
for date in dlist:
    yy=date[:4] ; mm=date[4:6] ; dd=date[6:8] ; hh=date[8:10]

    exp0grb=expsarch+'/'+explist[0]+'/'+date+'/pgbanl.gdas.'+date+'.grib2'
    exp1grb=expsarch+'/'+explist[1]+'/'+date+'/pgbanl.gdas.'+date+'.grib2'

    exp0=xa.open_dataset(exp0grb,filter_by_keys={"typeOfLevel": "isobaricInhPa"})
    exp1=xa.open_dataset(exp1grb,filter_by_keys={"typeOfLevel": "isobaricInhPa"})

    tmpexp0=exp0[expvar]
    tmpexp1=exp1[expvar]

    tmpexp0=tmpexp0.assign_coords(longitude=((tmpexp0.longitude+180)%360-180))
    tmpexp0=tmpexp0.sortby(tmpexp0.longitude)
    tmpexp1=tmpexp1.assign_coords(longitude=((tmpexp1.longitude+180)%360-180))
    tmpexp1=tmpexp1.sortby(tmpexp1.longitude)

    if (area!='Glb'):
       tmpexp0=tmpexp0.sel(latitude=slice(minlat,maxlat),longitude=slice(minlon,maxlon))
       tmpexp1=tmpexp1.sel(latitude=slice(minlat,maxlat),longitude=slice(minlon,maxlon))

    exp_tmpdiff=tmpexp1-tmpexp0

    if (date_cnt==0):
       exp_diff=exp_tmpdiff
    else:
       exp_diff=xa.concat((exp_diff,exp_tmpdiff),dim='time')

    date_cnt+=1

exp_diff_mean=exp_diff.mean(dim='time')

for pres in [850,500]:

    diff_da=exp_diff_mean.sel(isobaricInhPa=pres)

    melv_max=np.ceil( np.quantile(abs(diff_da),0.9) )
    melv_int=melv_max/10
    metcks=np.arange(-melv_max,melv_max+melv_int,melv_int)
    me_lvs=metcks
    clridx=[]
    for idx in np.linspace(2,254,metcks.size):
        clridx.append(int(idx))
    meclrmap=setup_cmap('BlueYellowRed',clridx)
    menorm = mpcrs.BoundaryNorm(me_lvs,len(clridx)+1,extend='both')
    
    fig,ax,gl=setupax_2dmap(cornerll,area,proj,lbsize=txsize)
    set_size(axe_w,axe_h,b=0.15,l=0.05,r=0.95,t=0.95)
    pltdata=diff_da
    cblabel='%s%s%s [%s]' %(expnlist[1],minussign,expnlist[0],units)
    cn=ax.contourf(pltdata.longitude,pltdata.latitude,pltdata,levels=me_lvs,
                   cmap=meclrmap,norm=menorm,extend='both')
    plt.colorbar(cn,orientation=cbori,fraction=cb_frac,
                 pad=cb_pad,ticks=me_lvs[::tkfreq],aspect=40,label=cblabel)
    outname='%s/%s_%s_%s_%s_Diff.%s_%s.%s' %(imgsavpath,expnlist[1],expnlist[0],pltvar,pres,sdate,edate,ffmt)
    if (fsave):
       print(outname)
       fig.savefig(outname,dpi=quality)
       plt.close()
