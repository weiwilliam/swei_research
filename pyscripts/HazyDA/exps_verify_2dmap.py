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
from utils import ndate,setup_cmap,find_cnlvs
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

outputpath=rootpath+'/vsERA5'
era5arch=rootarch+'/common/ERA5'
expsarch=rootarch+'/archive'

explist=['hazyda_ctrl','hazyda_aero']
expnlist=['CTL','AER']
pres_list=[850,700,500,200]

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

    exp0_tmperr=tmpexp0-tmpera
    exp1_tmperr=tmpexp1-tmpera

    if (date_cnt==0):
       exp0_err=exp0_tmperr
       exp1_err=exp1_tmperr
    else:
       exp0_err=xa.concat((exp0_err,exp0_tmperr),dim='time')
       exp1_err=xa.concat((exp1_err,exp1_tmperr),dim='time')

    date_cnt+=1

exp0_me=exp0_err.mean(dim='time')
exp1_me=exp1_err.mean(dim='time')

exp0_rmse=np.sqrt((exp0_err*exp0_err).mean(dim='time'))
exp1_rmse=np.sqrt((exp1_err*exp1_err).mean(dim='time'))

for pres in pres_list:

    me_da0=exp0_me.sel(isobaricInhPa=pres)
    me_da1=exp1_me.sel(isobaricInhPa=pres)
    me_da=xa.concat((me_da0,me_da1),dim='exps')
    rmse_da0=exp0_rmse.sel(isobaricInhPa=pres)
    rmse_da1=exp1_rmse.sel(isobaricInhPa=pres)
    rmse_da=xa.concat((rmse_da0,rmse_da1),dim='exps')

    me_lvs=find_cnlvs(me_da,eqside=1)
    clridx=[]
    for idx in np.linspace(2,254,me_lvs.size):
        clridx.append(int(idx))
    meclrmap=setup_cmap('BlueYellowRed',clridx)
    menorm = mpcrs.BoundaryNorm(me_lvs,len(clridx)+1,extend='both')

    dme_lvs=find_cnlvs((me_da1-me_da0),eqside=1)
    clridx=[]
    for idx in np.linspace(2,254,dme_lvs.size):
        clridx.append(int(idx))
    dmeclrmap=setup_cmap('BlueYellowRed',clridx)
    dmenorm = mpcrs.BoundaryNorm(dme_lvs,len(clridx)+1,extend='both')
    
    rmse_lvs=find_cnlvs(rmse_da)
    clridx=[]
    for idx in np.linspace(2,254,rmse_lvs.size):
        clridx.append(int(idx))
    rmseclrmap=setup_cmap('WhiteYellowOrangeRed',clridx)
    rmsenorm = mpcrs.BoundaryNorm(rmse_lvs,len(clridx)+1,extend='both')

    drmse_lvs=find_cnlvs((rmse_da1-rmse_da0),eqside=1)
    clridx=[]
    for idx in np.linspace(2,254,drmse_lvs.size):
        clridx.append(int(idx))
    drmseclrmap=setup_cmap('BlueYellowRed',clridx)
    drmsenorm = mpcrs.BoundaryNorm(drmse_lvs,len(clridx)+1,extend='both')

    fig,ax,gl=setupax_2dmap(cornerll,area,proj,lbsize=txsize)
    set_size(axe_w,axe_h,b=0.15,l=0.05,r=0.95)
    pltdata=exp0_me.sel(isobaricInhPa=pres)
    cblabel='%s Mean Error [%s]' %(expnlist[0],units)
    cn=ax.contourf(pltdata.longitude,pltdata.latitude,pltdata,levels=me_lvs,
                   cmap=meclrmap,norm=menorm,extend='both')
    plt.colorbar(cn,orientation=cbori,fraction=cb_frac,
                 pad=cb_pad,aspect=40,label=cblabel)
    outname='%s/%s_%s_%s_MeanError.%s_%s.%s' %(imgsavpath,expnlist[0],pltvar,pres,sdate,edate,ffmt)
    if (fsave):
       print(outname)
       fig.savefig(outname,dpi=quality)
       plt.close()

    fig,ax,gl=setupax_2dmap(cornerll,area,proj,lbsize=txsize)
    set_size(axe_w,axe_h,b=0.15,l=0.05,r=0.95)
    pltdata=exp1_me.sel(isobaricInhPa=pres)
    cblabel='%s Mean Error [%s]' %(expnlist[1],units)
    cn=ax.contourf(pltdata.longitude,pltdata.latitude,pltdata,levels=me_lvs,
                   cmap=meclrmap,norm=menorm,extend='both')
    plt.colorbar(cn,orientation=cbori,fraction=cb_frac,
                 pad=cb_pad,aspect=40,label=cblabel)
    outname='%s/%s_%s_%s_MeanError.%s_%s.%s' %(imgsavpath,expnlist[1],pltvar,pres,sdate,edate,ffmt)
    if (fsave):
       print(outname)
       fig.savefig(outname,dpi=quality)
       plt.close()

    fig,ax,gl=setupax_2dmap(cornerll,area,proj,lbsize=txsize)
    set_size(axe_w,axe_h,b=0.15,l=0.05,r=0.95)
    pltdata=exp1_me.sel(isobaricInhPa=pres)-exp0_me.sel(isobaricInhPa=pres)
    cblabel='%s%s%s Mean Error [%s]' %(expnlist[1],minussign,expnlist[0],units)
    cn=ax.contourf(pltdata.longitude,pltdata.latitude,pltdata,levels=dme_lvs,
                   cmap=dmeclrmap,norm=dmenorm,extend='both')
    plt.colorbar(cn,orientation=cbori,fraction=cb_frac,
                 pad=cb_pad,aspect=40,label=cblabel)
    title_str='(%.2f [%s])' %(pltdata.mean(),units)
    ax.set_title(title_str,loc='left')
    outname='%s/%s-%s_%s_%s_MeanError.%s_%s.%s' %(imgsavpath,expnlist[1],expnlist[0],pltvar,pres,sdate,edate,ffmt)
    if (fsave):
       print(outname)
       fig.savefig(outname,dpi=quality)
       plt.close()

    # RMSE
    fig,ax,gl=setupax_2dmap(cornerll,area,proj,lbsize=txsize)
    set_size(axe_w,axe_h,b=0.15,l=0.05,r=0.95)
    pltdata=exp0_rmse.sel(isobaricInhPa=pres)
    cblabel='%s Root-Mean-Square-Error [%s]' %(expnlist[0],units)
    cn=ax.contourf(pltdata.longitude,pltdata.latitude,pltdata,levels=rmse_lvs,
                   cmap=rmseclrmap,norm=rmsenorm,extend='both')
    plt.colorbar(cn,orientation=cbori,fraction=cb_frac,
                 pad=cb_pad,aspect=40,label=cblabel)
    outname='%s/%s_%s_%s_RMSE.%s_%s.%s' %(imgsavpath,expnlist[0],pltvar,pres,sdate,edate,ffmt)
    if (fsave):
       print(outname)
       fig.savefig(outname,dpi=quality)
       plt.close()

    fig,ax,gl=setupax_2dmap(cornerll,area,proj,lbsize=txsize)
    set_size(axe_w,axe_h,b=0.15,l=0.05,r=0.95)
    pltdata=exp1_rmse.sel(isobaricInhPa=pres)
    cblabel='%s Root-Mean-Square-Error [%s]' %(expnlist[1],units)
    cn=ax.contourf(pltdata.longitude,pltdata.latitude,pltdata,levels=rmse_lvs,
                   cmap=rmseclrmap,norm=rmsenorm,extend='both')
    plt.colorbar(cn,orientation=cbori,fraction=cb_frac,
                 pad=cb_pad,aspect=40,label=cblabel)
    outname='%s/%s_%s_%s_RMSE.%s_%s.%s' %(imgsavpath,expnlist[1],pltvar,pres,sdate,edate,ffmt)
    if (fsave):
       print(outname)
       fig.savefig(outname,dpi=quality)
       plt.close()

    fig,ax,gl=setupax_2dmap(cornerll,area,proj,lbsize=txsize)
    set_size(axe_w,axe_h,b=0.15,l=0.05,r=0.95)
    pltdata=exp1_rmse.sel(isobaricInhPa=pres)-exp0_rmse.sel(isobaricInhPa=pres)
    cblabel='%s%s%s Root-Mean-Square-Error [%s]' %(expnlist[1],minussign,expnlist[0],units)
    cn=ax.contourf(pltdata.longitude,pltdata.latitude,pltdata,levels=drmse_lvs,
                   cmap=drmseclrmap,norm=drmsenorm,extend='both')
    plt.colorbar(cn,orientation=cbori,fraction=cb_frac,
                 pad=cb_pad,aspect=40,label=cblabel)
    title_str='(%.2f [%s])' %(pltdata.mean(),units)
    ax.set_title(title_str,loc='left')
    outname='%s/%s-%s_%s_%s_RMSE.%s_%s.%s' %(imgsavpath,expnlist[1],expnlist[0],pltvar,pres,sdate,edate,ffmt)
    if (fsave):
       print(outname)
       fig.savefig(outname,dpi=quality)
       plt.close()
