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
elif (machine=='Cheyenne'):
    rootpath='/glade/work/swei/output/images'
    rootarch='/scratch2/BMC/gsd-fv3-dev/Shih-wei.Wei/ResearchData'
    rootgit='/glade/u/home/swei/research'
elif (machine=='S4'):
    rootpath='/data/users/swei/Images'
    rootarch='/scratch/users/swei/ncdiag'
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
txsize=12
mpl.rc('axes',titlesize=txsize,labelsize=txsize)
mpl.rc('xtick',labelsize=txsize)
mpl.rc('ytick',labelsize=txsize)
mpl.rc('legend',fontsize='x-large')
axe_w=6; axe_h=3; quality=300

if (machine=='Cheyenne'):
   inputpath='/glade/work/dfgrogan/UFS/WM_DTAER/AER'
elif (machine=='S4'):
   inputpath='/data/users/swei/common/MERRA2'

m2tag='inst3_3d_aer_Nv'
tkfreq=2

sdate=2020081300
edate=2020092118
hint=6
pltvar='total'
area='Glb'
pltall=0 # 0: total only 1: sub species included
pltsnap=0
pltave=1
plt_pts=0
#
if (plt_pts):
   pts_lat=[17.98]
   pts_lon=[-60.72]
   ptsize=20
# 
clridx=[0,11,20,29,38,47,56,65,74,83,92,101,110,119,128]
cn_cmap=setup_cmap('MPL_YlOrBr',clridx)
cnlvs=[0      ,  2.5e-5,   5e-5, 7.5e-5,
       1.25e-4,  1.5e-4,1.75e-4,   2e-4,
          3e-4,    4e-4,   5e-4,   6e-4,
          8e-4,    1e-3,   3e-3,   5e-3]
#cnlvs=[0      ,  2.5e-6,   5e-6, 7.5e-6,
#       1.25e-5,  1.5e-5,1.75e-5,   2e-5,
#          3e-5,    4e-5,   5e-5,   6e-5,
#          8e-5,    1e-4,   3e-4,   5e-4]
outputpath=rootpath+'/Dataset/MERRA-2/2dmap/'+area
if ( not os.path.exists(outputpath) ):
    os.makedirs(outputpath)
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

cblb='Column mass density [$\\times\mathrm{10}^{-%i}$ $\mathrm{kg\cdot m^{-2}}$]' %(int(np.log10(scalef)))
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

    print('Processing '+infile,flush=1)
   
    ds=xa.open_dataset(infile)
    if (multi_time):
       ds=ds.sel(time=dates[dates_count])
    else:
       ds=ds.sel(time=ds.time[0])
   
    delp=ds.DELP
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

    cmass=conc.sum(dim='lev')
    
    if (not pltall):
       nplotlist=[nvars]
    else:
       nplotlist=np.arange(nvars+1)
   
    if (pltsnap): 
       for n in nplotlist:
           if (n<nvars):
              title='%s column mass density' %(varlst[n])
              outname='%s/%s_%s_cmass.%s.png'  %(outputpath,area,varlst[n],date)
           else:
              title='%s column mass density' %(varname)
              outname='%s/%s_%s_all_cmass.%s.png' %(outputpath,area,pltvar,date)
       
           pltdata=cmass[n,:,:]
           fig,ax,gl=setupax_2dmap(cornerll,area,proj,lbsize=txsize)
           set_size(axe_w,axe_h,b=0.15,l=0.05,r=0.95,t=0.95)
           cn=ax.contourf(pltdata.lon,pltdata.lat,pltdata*scalef,levels=cnlvsarr,cmap=cn_cmap,norm=norm)
           #ax.set_title(title)
           plt.colorbar(cn,ax=ax,orientation='horizontal',
                        format='%.2f',fraction=0.04,aspect=40,pad=0.08,label=cblb)
           if (plt_pts):
              sc=ax.scatter(pts_lon,pts_lat,s=ptsize,c='w',marker='x')
           print(outname)
           fig.savefig(outname,dpi=quality)
           plt.close()
 
    if (pltave):
       if (dates_count==0):
          cmass_all=cmass
       else:
          cmass_all=xa.concat((cmass_all,cmass),dim='time')

    dates_count+=1

if (pltave):
   cmass_mean=cmass_all.mean(dim='time')
   for n in nplotlist:
       if (n<nvars):
          title='%s column mass density' %(varlst[n])
          outname='%s/%s_%s_cmass.%s_%s.png'  %(outputpath,area,varlst[n],sdate,edate)
       else:
          title='%s column mass density' %(varname)
          outname='%s/%s_%s_all_cmass.%s_%s.png' %(outputpath,area,pltvar,sdate,edate)

       pltdata=cmass_mean.sel(bins=n)
       fig,ax,gl=setupax_2dmap(cornerll,area,proj,lbsize=txsize)
       set_size(axe_w,axe_h,b=0.15,l=0.05,r=0.95,t=0.95)
       cn=ax.contourf(pltdata.lon,pltdata.lat,pltdata*scalef,levels=cnlvsarr,cmap=cn_cmap,norm=norm)
       #ax.set_title(title)
       plt.colorbar(cn,ax=ax,orientation='horizontal',
                    format='%.2f',fraction=0.04,aspect=40,pad=0.08,label=cblb)
       if (plt_pts):
          sc=ax.scatter(pts_lon,pts_lat,s=ptsize,c='w',marker='x')
       print(outname)
       fig.savefig(outname,dpi=quality)
       plt.close()
