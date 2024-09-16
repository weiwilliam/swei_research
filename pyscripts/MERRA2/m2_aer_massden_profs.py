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
        machine='Cheyenne'
    elif (os.path.exists('/cardinal')):
        rootpath='/data/users/swei/Images'
        rootarch='/scratch2/BMC/gsd-fv3-dev/Shih-wei.Wei/ResearchData'
        rootgit='/home/swei/research'
        machine='S4'
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
txsize=12
mpl.rc('axes',titlesize=txsize,labelsize=txsize)
mpl.rc('xtick',labelsize=txsize)
mpl.rc('ytick',labelsize=txsize)
mpl.rc('legend',fontsize='large')
fsave=1; figfmt='png'; quality=300
axe_w=3; axe_h=5

if (machine=='Cheyenne'):
   inputpath='/glade/work/dfgrogan/UFS/WM_DTAER/AER'
elif (machine=='S4'):
   inputpath='/data/users/swei/common/MERRA2'

outputpath=rootpath+'/Dataset/MERRA-2/area_profs'
if ( not os.path.exists(outputpath) ):
    os.makedirs(outputpath)

sdate=2020091112
edate=2020091112
hint=6
pltvar='carbon'
area='Y20Smk1'
pltall=0 # 0: total only 1: sub species included
m2tag='inst3_3d_aer_Nv'
tkfreq=2

# 
# Constant configuration
proj=ccrs.PlateCarree()
grav=9.80665e+0

# Setup plotting area
minlon, maxlon, minlat, maxlat, crosszero, cyclic=setarea.setarea(area)
print(minlat,maxlat,minlon,maxlon,crosszero,cyclic)
if (area=='Glb'):
   minlon=-180. ; maxlon=180.
else:
   minlon=(minlon+180)%360-180
   maxlon=(maxlon+180)%360-180
cornerll=[minlat,maxlat,minlon,maxlon]

if (pltvar=='dust'):
   varlst=['DU001','DU002','DU003','DU004','DU005']
   varname='dust'
   scalef=1e4
   cnscale=1e4
elif (pltvar=='seas'):
   varlst=['SS001','SS002','SS003','SS004','SS005']
   varname='sea salt'
   scalef=1e4
   cnscale=1e4
elif (pltvar=='carbon'):
   varlst=['OCPHILIC','OCPHOBIC','BCPHILIC','BCPHOBIC'] 
   varname='carbonaceous'
   scalef=1e5
   cnscale=1e4
elif (pltvar=='sulf'):
   varlst=['SO4']
   varname='sulfate'
   scalef=1e5
   cnscale=1e4
elif (pltvar=='total'):
   varlst=['DU001','DU002','DU003','DU004','DU005',
           'SS001','SS002','SS003','SS004','SS005',
           'OCPHILIC','OCPHOBIC','BCPHILIC','BCPHOBIC',
           'SO4']
   varname='total'
   scalef=1e4
   cnscale=1e4

#cblb='Column mass density [%.0e $\mathrm{kg\cdot m^{-2}}$]' %(1/scalef)
#cnlvsarr=np.array(cnlvs)*cnscale
#norm = mpcrs.BoundaryNorm(cnlvsarr,len(cnlvsarr))

nvars=len(varlst)

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

dates_count=0
for date in dlist:
    yy=date[:4] ; mm=date[4:6] ; dd=date[6:8] ; hh=date[8:10]
    pdy=date[:8]
    if (yy=='2020' and mm=='09'):
       m2ind='401'
    else:
       m2ind='400'
# MERRA2_401.inst3_3d_aer_Nv.20200916_12Z.nc4
    if (machine=='Cheyenne'):
       infile=inputpath+'/MERRA2_'+m2ind+'.'+m2tag+'.'+pdy+'_'+hh+'Z.nc4'
       multi_time=0
    if (machine=='S4'):
       infile=inputpath+'/'+yy+'/'+mm+'/MERRA2_'+m2ind+'.'+m2tag+'.'+pdy+'.nc4'
       multi_time=1
   
    ds=xa.open_dataset(infile)
    if (multi_time):
       ds=ds.sel(time=dates[dates_count])
    else:
       ds=ds.sel(time=ds.time[0])
    ds=ds.sel(lon=slice(minlon,maxlon),lat=slice(minlat,maxlat))

    delp=ds.DELP
    kgkg_kgm2=delp/grav
    pres=xa.zeros_like(delp)
    for i in np.arange(ds.lev.size):
        pres[i,:,:]=ds.PS-delp[i:,:,:].sum(dim='lev')
    pres=pres/100.

    for var in varlst:
        tmp=ds[var]*kgkg_kgm2
        if (var==varlst[0]):
           conc=tmp
           total=tmp
        else:
           conc=xa.concat((conc,tmp),dim='bins')
           total=total+tmp

    conc=xa.concat((conc,total),dim='bins')

# Get the column integral
    cmass=np.sum(conc,axis=1)
    
    if (not pltall):
       nplotlist=[nvars]
    else:
       nplotlist=np.arange(nvars+1)
#    
    xlb='Mass Density [$\mathrm{kg\cdot m^{-2}}$]'
    ylb='Pressure [hPa]'
    for n in nplotlist:
        third_quantile=cmass[n].quantile(0.75)
        if (n<nvars):
           outname='%s/%s_%s_prof.%s.%s'  %(outputpath,area,varlst[n],date,figfmt)
        else:
           outname='%s/%s_%s_all_prof.%s.%s' %(outputpath,area,pltvar,date,figfmt)

        profdata=conc.sel(bins=n)
        msk=xa.ones_like(profdata,dtype='bool')
        pltmsk=(msk)&(cmass[n]>third_quantile)
        area_data=xa.where(msk,profdata,np.nan)
        area_pres=xa.where(msk,pres,np.nan)
        pltprof=area_data.mean(axis=(1,2),skipna=1)
        pltpres=area_pres.mean(axis=(1,2),skipna=1)
        pltprofs=np.reshape(area_data.values,(area_data.lev.size,area_data.lon.size*area_data.lat.size))

        fig,ax=plt.subplots()
        set_size(axe_w,axe_h,l=0.2)
        ax.plot(pltprofs,pltpres,'-',c='tab:red',linewidth=0.5,alpha=0.6)
        ax.plot(pltprof,pltpres,'-',c='k',linewidth=1.5)
        ax.set_yscale('log')
        ax.set_xlabel(xlb)
        ax.set_ylabel(ylb)
        ax.set_ylim(200,1025)
        ax.get_yaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
        ax.get_yaxis().set_minor_formatter(mpl.ticker.ScalarFormatter())
        ax.invert_yaxis()
        #ax.legend(lglst)
        ax.grid(which='both')
#    
        if (fsave):
           print(outname)
           fig.savefig(outname,dpi=quality)
           plt.close()
 
    dates_count+=1
