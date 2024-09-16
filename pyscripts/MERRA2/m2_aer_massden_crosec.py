import sys, os, platform
machine='S4'
if (machine=='MBP'):
    rootpath='/Users/weiwilliam'
    rootarch='/Volumes/WD2TB/ResearchData'
elif (machine=='Desktop'):
    rootpath='F:\GoogleDrive_NCU\Albany'
    rootarch='F:\ResearchData'
    rootgit='F:\GitHub\swei_research'
elif (machine=='Hera'):
    rootpath='/scratch2/BMC/gsd-fv3-dev/Shih-wei.Wei'
    rootarch='/scratch2/BMC/gsd-fv3-dev/Shih-wei.Wei/ResearchData'
    rootgit='/home/Shih-wei.Wei/research'
elif (machine=='Cheyenne')):
    rootpath='/glade/work/swei/output/images'
    rootarch='/scratch2/BMC/gsd-fv3-dev/Shih-wei.Wei/ResearchData'
    rootgit='/glade/u/home/swei/research'
elif (machine='S4')):
    rootpath='/data/users/swei/Images'
    rootarch='/scratch/users/swei/ncdiag'
    rootgit='/home/swei/research'
sys.path.append(rootgit+'/pyscripts/functions')
from utils import setup_cmap, ndate, lat_ns, lon_we
from plot_utils import setupax_2dmap,set_size
import setuparea as setarea
import xarray as xa
import numpy as np
import pandas as pd
from datetime import datetime, timedelta
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mpcrs
from matplotlib import ticker
import cartopy.crs as ccrs

# Plotting setup
mpl.rc('axes',titlesize=18,labelsize=16)
mpl.rc('xtick',labelsize=16)
mpl.rc('ytick',labelsize=16)
mpl.rc('legend',fontsize='xx-large')
axe_w=10; axe_h=5
quality=300
deg_sym=u'\u00B0'

if (machine=='Cheyenne'):
   inputpath='/glade/work/dfgrogan/UFS/WM_DTAER/AER'
elif (machine=='S4'):
   inputpath='/data/users/swei/common/MERRA2'

outputpath=rootpath+'/Dataset/MERRA-2/Crossec'
if ( not os.path.exists(outputpath) ):
    os.makedirs(outputpath)

sdate=2020062212
edate=2020062212
hint=6
pltvar='dust'
area='r2o1'
pltall=0 # 0: total only 1: sub species included
m2tag='inst3_3d_aer_Nv'
tkfreq=2
crotype='meridional' # zonal/meridional

# 
clridx=[0,11,20,29,38,47,56,65,74,83,92,101,110,119,128]
cn_cmap=setup_cmap('MPL_YlOrBr',clridx)
#cnlvs=[0      ,  2.5e-5,   5e-5, 7.5e-5,
#       1.25e-4,  1.5e-4,1.75e-4,   2e-4,
#          3e-4,    4e-4,   5e-4,   6e-4,
#          8e-4,    1e-3,   3e-3,   5e-3]
cnlvs=[0      ,  2.5e-6,   5e-6, 7.5e-6,
       1.25e-5,  1.5e-5,1.75e-5,   2e-5,
          3e-5,    4e-5,   5e-5,   6e-5,
          8e-5,    1e-4,   3e-4,   5e-4]
prescnlvs=['500','600','700','800','900','1000']

# Constant configuration
proj=ccrs.PlateCarree()
grav=9.80665e+0

# Setup plotting area
minlon, maxlon, minlat, maxlat, crosszero, cyclic=setarea.setarea(area)
print(minlat,maxlat,minlon,maxlon,crosszero,cyclic)
if (area=='Glb'):
   minlon=-180. ; maxlon=180.
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

cblb='Mass density [%.0e $\mathrm{kg\cdot m^{-2}}$]' %(1/scalef)
cnlvsarr=np.array(cnlvs)*cnscale
norm = mpcrs.BoundaryNorm(cnlvsarr,len(cnlvsarr))

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

    if (area!='Glb'):
       ds=ds.sel(lon=slice(minlon,maxlon),lat=slice(minlat,maxlat))
   
    delp=ds.DELP
    pres=xa.zeros_like(delp)
    for i in np.arange(ds.lev.size):
        pres[i,:,:]=ds.PS-delp[i:,:,:].sum(dim='lev')
    pres=pres/100.
    domain_ave_pres=pres.mean(dim=('lon','lat'))
    
    kgkg_kgm2=delp/grav

    for var in varlst:
        tmp=ds[var]*kgkg_kgm2

        if (var==varlst[0]):
           conc=tmp
           total=tmp
        else:
           conc=xa.concat((conc,tmp),dim='bins')
           total=total+tmp

    conc=xa.concat((conc,total),dim='bins')

# Get the zonal or meridional mean
    if (crotype=='zonal'):
       cross=conc.mean(dim='lon')
       cross_pres=pres.mean(dim='lon')
    elif (crotype=='meridional'):
       cross=conc.mean(dim='lat')
       cross_pres=pres.mean(dim='lat')
    
    if (not pltall):
       nplotlist=[nvars]
    else:
       nplotlist=np.arange(nvars+1)
    
    for n in nplotlist:
        if (n<nvars):
           title='%s mass density' %(varlst[n])
           outname='%s/%s_%s_%s.%s.png'  %(outputpath,area,varlst[n],crotype,date)
        else:
           title='%s mass density' %(varname)
           outname='%s/%s_%s_all_%s.%s.png' %(outputpath,area,pltvar,crotype,date)
    
        pltdata=cross[n,:,:]
        fig,ax=plt.subplots()
        set_size(axe_w,axe_h,b=0.13,l=0.10,r=0.95,t=0.95)
        if (crotype=='zonal'):
           cn=ax.contourf(pltdata.lat,pltdata.lev,pltdata*scalef,levels=cnlvsarr,cmap=cn_cmap,norm=norm)
           prescn=ax.contour(pltdata.lat,pltdata.lev,cross_pres,
                             linestyles='dashed',levels=prescnlvs,colors='k')
           x_format=ticker.FuncFormatter(lambda x, pos: "%s" %(lat_ns(x))) 
        elif (crotype=='meridional'):
           cn=ax.contourf(pltdata.lon,pltdata.lev,pltdata*scalef,levels=cnlvsarr,cmap=cn_cmap,norm=norm)
           prescn=ax.contour(pltdata.lon,pltdata.lev,cross_pres,
                             linestyles='dashed',levels=prescnlvs,colors='k')
           x_format=ticker.FuncFormatter(lambda x, pos: "%s" %(lon_we(x)))
        ax.clabel(prescn, prescn.levels, inline=True, fmt='%.0f hPa', fontsize=10)
        ax.set_title(title)
        ax.xaxis.set_major_formatter(x_format)
        ax.invert_yaxis()
        ax.set_ylabel('Model level')
        plt.colorbar(cn,ax=ax,orientation='horizontal',ticks=cnlvsarr[::tkfreq],
                     format='%.2f',fraction=0.045,aspect=40,pad=0.08,label=cblb)
        print(outname)
        fig.savefig(outname,dpi=quality)
        plt.close()
 
    dates_count+=1
