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
        machine='hera'
        rootpath='/scratch2/BMC/gsd-fv3-dev/Shih-wei.Wei'
        rootarch='/scratch2/BMC/gsd-fv3-dev/Shih-wei.Wei/ResearchData'
        rootgit='/home/Shih-wei.Wei/research'
    elif (os.path.exists('/glade')):
        machine='cheyenne'
        rootpath='/glade/work/swei/output/images'
        rootarch='/scratch2/BMC/gsd-fv3-dev/Shih-wei.Wei/ResearchData'
        rootgit='/glade/u/home/swei/research'
    elif (os.path.exists('/scratch')):
        machine='s4'
        rootpath='/data/users/swei/Images'
        rootarch='/scratch2/BMC/gsd-fv3-dev/Shih-wei.Wei/ResearchData'
        rootgit='/home/swei/research'
print('%s %s' % (os_name,machine))
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
axe_w=10; axe_h=5
quality=300
tkfreq=1

inputpath='/data/users/swei/archive'
#inputpath='/scratch/users/swei/comrot'
outputpath=rootpath+'/Ens_size/Timeseries'
if ( not os.path.exists(outputpath) ):
    os.makedirs(outputpath)

sdate=2020060100
edate=2020061218
hint=6
pltvar='total'
area='Glb'
explist=['ctrl','10ens','30ens']
evallist=['pressfc','dpres','tmp','spfh','ugrd','vgrd']

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

for i in np.arange(tnum-1):
    print(dlist[i])
    eidx=0
    for exp in explist:
        if (exp=='ctrl'):
           numens=20
        else:
           numens=int(exp[:2])

        adate=dlist[i+1]
        print(adate)
        a_yy=adate[:4] ; a_mm=adate[4:6] ; a_dd=adate[6:8] ; a_hh=adate[8:10]
        a_pdy=adate[:8]
        ensres_anl=inputpath+'/hazyda_'+exp+'/'+adate+'/gdas.t'+a_hh+'z.atmanl.ensres.nc'
        #ensres_anl=inputpath+'/hazyda_'+exp+'/gdas.'+a_pdy+'/'+a_hh+'/atmos/gdas.t'+a_hh+'z.atmanl.ensres.nc'
        anl=xa.open_dataset(ensres_anl)
        ensanl=anl.get(evallist)

        gdate=dlist[i]
        g_yy=gdate[:4] ; g_mm=gdate[4:6] ; g_dd=gdate[6:8] ; g_hh=gdate[8:10]
        g_pdy=gdate[:8]

        for nens in np.arange(1,numens+1):
            nens_str='%.3i' %(nens)
            ens6hr=inputpath+'/hazyda_'+exp+'/'+gdate+'/enkf/mem'+nens_str+'/gdas.t'+g_hh+'z.atmf006.nc'
            #ens6hr=inputpath+'/hazyda_'+exp+'/enkfgdas.'+g_pdy+'/'+g_hh+'/atmos/mem'+nens_str+'/gdas.t'+g_hh+'z.atmf006.nc'
            ensmem=xa.open_dataset(ens6hr)
            ensmem=ensmem.get(evallist)
            mean_err2=(ensmem-ensanl)*(ensmem-ensanl)/numens
            if (nens==1):
               sum_mean_err2=mean_err2
            else:
               sum_mean_err2+=mean_err2 
        rmse_exp=np.sqrt(sum_mean_err2)

        enssprd=inputpath+'/hazyda_'+exp+'/'+gdate+'/enkf/gdas.t'+g_hh+'z.atmf006.ensspread.nc'
        #enssprd=inputpath+'/hazyda_'+exp+'/enkfgdas.'+g_pdy+'/'+g_hh+'/atmos/gdas.t'+g_hh+'z.atmf006.ensspread.nc'
        sprd_exp=xa.open_dataset(enssprd)
        sprd_exp=sprd_exp.get(evallist)

        #ratio_tmp=rmse/sprd
        #ratio_tmp=sprd/rmse
        #if (i==0):
        #   rmse_exp=rmse
        #   sprd_exp=sprd
        #else:
        #   rmse_exp=xa.concat((rmse_exp,rmse),dim='time')
        #   sprd_exp=xa.concat((sprd_exp,sprd),dim='time')

        rmse_exp=rmse_exp.assign_coords({'exps':exp})
        sprd_exp=sprd_exp.assign_coords({'exps':exp})
        if (eidx==0):
           rmse=rmse_exp
           sprd=sprd_exp
        else:
           rmse=xa.concat((rmse,rmse_exp),dim='exps')
           sprd=xa.concat((sprd,sprd_exp),dim='exps')
        eidx+=1

    inputpath='/data/users/swei/archive'
    outputpath=inputpath+'/Ens_size'
    if ( not os.path.exists(outputpath) ):
        os.makedirs(outputpath)
    
    rmse_fname=outputpath+'/rmse.'+adate+'.nc4'
    ratio.to_netcdf(rmse_fname)
    sprd_fname=outputpath+'/sprd.'+adate+'.nc4'
    ratio.to_netcdf(sprd_fname)
