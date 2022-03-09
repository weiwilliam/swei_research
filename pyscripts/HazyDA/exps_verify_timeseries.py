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
    rootpath='/data/users/swei'
    rootarch='/scratch/users/swei/ncdiag'
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
import scipy.stats

tlsize=12 ; lbsize=10
mpl.rc('axes', titlesize=tlsize,labelsize=lbsize)
mpl.rc('xtick',labelsize=lbsize)
mpl.rc('ytick',labelsize=lbsize)
mpl.rc('legend',fontsize='small')
fsave=1 ; ffmt='png' ; ptsize=4
axe_w=3 ; axe_h=3 ; quality=300
minussign=u'\u2212'

outputpath=rootpath+'/vsERA5'
era5arch=rootarch+'/common/ERA5'
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

imgsavpath=outputpath+'/'+area+'/2dmap'
if ( not os.path.exists(imgsavpath) ):
   os.makedirs(imgsavpath)

syy=int(str(sdate)[:4]); smm=int(str(sdate)[4:6])
sdd=int(str(sdate)[6:8]); shh=int(str(sdate)[8:10])
eyy=int(str(edate)[:4]); emm=int(str(edate)[4:6])
edd=int(str(edate)[6:8]); ehh=int(str(edate)[8:10])

date1 = datetime(syy,smm,sdd,shh)
date2 = datetime(eyy,emm,edd,ehh)
delta = timedelta(hours=hint)
dates = pd.date_range(start=date1, end=date2, freq=delta)

xdate2= date2+delta
xdates= mdates.drange(date1, xdate2, delta)

rule = rrulewrapper(DAILY, byhour=6, interval=5)
loc = RRuleLocator(rule)
formatter = DateFormatter('%Y%h %n %d %Hz')

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
    eragrb=era5arch+'/'+yy+'/'+mm+'/era5_'+date+'.grib'

    exp0grb=expsarch+'/'+explist[0]+'/'+date+'/pgbanl.gdas.'+date+'.grib2'
    exp1grb=expsarch+'/'+explist[1]+'/'+date+'/pgbanl.gdas.'+date+'.grib2'

    if (os.path.exists(eragrb)):
       era=xa.open_dataset(eragrb)
    
    exp0=xa.open_dataset(exp0grb,filter_by_keys={"typeOfLevel": "isobaricInhPa"})
    exp1=xa.open_dataset(exp1grb,filter_by_keys={"typeOfLevel": "isobaricInhPa"})

    tmpera=era[pltvar]
    if (pltvar=='z'):
       tmpera=tmpera/grav
    tmpexp0=exp0[expvar]
    tmpexp1=exp1[expvar]
    tmpera=tmpera.interp_like(tmpexp0)

    tmpera=tmpera.assign_coords(longitude=((tmpera.longitude+180)%360-180))
    tmpera=tmpera.sortby(tmpera.longitude)
    tmpexp0=tmpexp0.assign_coords(longitude=((tmpexp0.longitude+180)%360-180))
    tmpexp0=tmpexp0.sortby(tmpexp0.longitude)
    tmpexp1=tmpexp1.assign_coords(longitude=((tmpexp1.longitude+180)%360-180))
    tmpexp1=tmpexp1.sortby(tmpexp1.longitude)

    if (area!='Glb'):
       tmpera=tmpera.sel(latitude=slice(minlat,maxlat),longitude=slice(minlon,maxlon))
       tmpexp0=tmpexp0.sel(latitude=slice(minlat,maxlat),longitude=slice(minlon,maxlon))
       tmpexp1=tmpexp1.sel(latitude=slice(minlat,maxlat),longitude=slice(minlon,maxlon))

    lat_weight=np.cos(np.deg2rad(tmpera.latitude))
    exp0_err=tmpexp0-tmpera
    exp1_err=tmpexp1-tmpera

    exp0_me=(exp0_err.weighted(lat_weight)).mean(('longitude','latitude'))
    exp1_me=(exp1_err.weighted(lat_weight)).mean(('longitude','latitude'))

    me0_da=xa.Dataset({pltvar:(['preslv'],exp0_me.data)},
                      coords={'preslv':exp0_me.isobaricInhPa.data})
    me1_da=xa.Dataset({pltvar:(['preslv'],exp1_me.data)},
                      coords={'preslv':exp1_me.isobaricInhPa.data})

    exp0_sqerr_w=(exp0_err*exp0_err).weighted(lat_weight)
    exp1_sqerr_w=(exp1_err*exp1_err).weighted(lat_weight)

    exp0_rmse=np.sqrt( (exp0_sqerr_w).mean(('longitude','latitude')) )
    exp1_rmse=np.sqrt( (exp1_sqerr_w).mean(('longitude','latitude')) )

    rmse0_da=xa.Dataset({pltvar:(['preslv'],exp0_rmse.data)},
                        coords={'preslv':exp0_rmse.isobaricInhPa.data})
    rmse1_da=xa.Dataset({pltvar:(['preslv'],exp1_rmse.data)},
                        coords={'preslv':exp1_rmse.isobaricInhPa.data})

    if (date_cnt==0):
       me0_ts_da=me0_da
       me1_ts_da=me1_da
       rmse0_ts_da=rmse0_da
       rmse1_ts_da=rmse1_da
    else:
       me0_ts_da=xa.concat((me0_ts_da,me0_da),dim='time')
       me1_ts_da=xa.concat((me1_ts_da,me1_da),dim='time')
       rmse0_ts_da=xa.concat((rmse0_ts_da,rmse0_da),dim='time')
       rmse1_ts_da=xa.concat((rmse1_ts_da,rmse1_da),dim='time')

    date_cnt+=1


plt_me=(xa.concat((me0_ts_da,me1_ts_da),dim='exps'))
plt_rmse=(xa.concat((rmse0_ts_da,rmse1_ts_da),dim='exps'))
preslv=plt_me.preslv.data
biasprof=plt_me.mean(dim='time')
rmseprof=plt_rmse.mean(dim='time')
diffprof=xa.concat((biasprof.diff(dim='exps'),rmseprof.diff(dim='exps')),dim='type')
diffprof=diffprof.sel(exps=0)

fig,ax=plt.subplots(1,3,sharey=True)
fig.subplots_adjust(wspace=0.05)
ax[0].set_prop_cycle(color=['red','blue'])
ax[1].set_prop_cycle(color=['red','blue'])
ax[2].set_prop_cycle(color=['k','k'],linestyle=['-','--'])
ax[0].plot(biasprof[pltvar].data.swapaxes(0,1),preslv,'-')
ax[0].set_title('BIAS [%s]'%(units))
ax[1].plot(rmseprof[pltvar].data.swapaxes(0,1),preslv,'--')
ax[1].set_title('RMSE [%s]'%(units))
ax[2].plot(diffprof[pltvar].data.swapaxes(0,1),preslv)
ax[2].set_title(expnlist[1]+minussign+expnlist[0])
ax[0].invert_yaxis()
fig.suptitle('%s [%s]' %(pltvar.upper(),units))
ax[0].legend(expnlist)
ax[2].legend(['BIAS','RMSE'])
ax[0].grid()
ax[1].grid()
ax[2].grid()
if (fsave):
    filename=imgsavpath+'/%s_%s_%s_%s_BIASRMS.png' %(area,pltvar,expnlist[0],expnlist[1])
    fig.savefig(filename,dpi=quality)
    plt.close()

for pres in [850,500]:
    fig,ax=plt.subplots(2,1,sharex=True,figsize=(9,3.8))
    fig.subplots_adjust(hspace=0.15)
    for a in np.arange(2):
        ax[a].set_prop_cycle(color=['red','blue'])
        ax[a].grid()
    
    plt_lv_me=plt_me[pltvar].sel(preslv=pres).data.swapaxes(0,1)
    plt_lv_rmse=plt_rmse[pltvar].sel(preslv=pres).data.swapaxes(0,1)
    ax[0].xaxis.set_major_locator(loc)
    ax[0].xaxis.set_major_formatter(formatter)
    ax[0].xaxis.set_tick_params(rotation=30, labelsize=10)
    ax[0].set_title('%s [%s] %.1f [hPa]' %(pltvar.upper(),units,pres),loc='left')
    ax[0].plot_date(dates,plt_lv_rmse,'--')
    ax[0].set_ylabel('RMSE [%s]'%(units))
    ax[1].plot_date(dates,plt_lv_me,'-')
    ax[1].set_ylabel('BIAS [%s]'%(units))
    lglist=np.zeros((2,2),dtype='<U30')
    for ex in np.arange(2):
        lglist[0,ex]=expnlist[ex]+'(%.2f)' %(plt_lv_rmse[:,ex].mean())
        lglist[1,ex]=expnlist[ex]+'(%.2f)' %(plt_lv_me[:,ex].mean())
    ax[0].legend(lglist[0,:])
    ax[1].legend(lglist[1,:])

    if (fsave):
        fig.savefig(imgsavpath+'/%s_%s_%s_%s_%i_BIASRMS.png'
                    %(area,pltvar,expnlist[0],expnlist[1],pres), dpi=200)
        plt.close()


