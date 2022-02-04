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
from pydap.cas.urs import setup_session

# Plotting setup
mpl.rc('axes',titlesize=18,labelsize=16)
mpl.rc('xtick',labelsize=16)
mpl.rc('ytick',labelsize=16)
mpl.rc('legend',fontsize='xx-large')
axe_w=10; axe_h=5
quality=300

#inputpath='/glade/work/dfgrogan/UFS/WM_DTAER/AER'
m2_url = 'https://goldsmr4.gesdisc.eosdis.nasa.gov/opendap/MERRA2/M2T1NXAER.5.12.4'
#/2020/09/MERRA2_401.tavg1_2d_aer_Nx.20200901.nc4'
outputpath=rootpath+'/Dataset/MERRA-2/2dmap/AOD'
if ( not os.path.exists(outputpath) ):
    os.makedirs(outputpath)

sdate=2020082200
edate=2020082218
hint=6
pltvar='DUEXTTAU'
area='Glb'
pltall=0 # 0: total only 1: sub species included
#m2tag='inst3_2d_gas_Nx' #Total AOD
m2tag='tavg1_2d_aer_Nx' # Speciated AOD
tkfreq=1

#
cnlvs=np.array((0., 0.05, 0.1, 0.15, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 2.5))
clridx=np.array((0,2,4,7,8,9,10,12,14,15,16,17,18))
clrmap=setup_cmap('precip3_16lev',clridx)
aer_norm = mpcrs.BoundaryNorm(cnlvs,len(clridx),extend='both')
cblb=pltvar

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

p_pdy='00000000'
dates_count=0
for date in dlist:
    yy=date[:4] ; mm=date[4:6] ; dd=date[6:8] ; hh=date[8:10]
    pdy=date[:8]
    if (yy=='2020' and mm=='09'):
       m2ind='401'
    else:
       m2ind='400'
#/2020/09/MERRA2_401.tavg1_2d_aer_Nx.20200901.nc4'
    if (pdy!=p_pdy):
       in_url=m2_url+'/'+yy+'/'+mm+'/MERRA2_'+m2ind+'.'+m2tag+'.'+pdy+'.nc4'
       session = setup_session('username', 'password', check_url=in_url)
       store = xa.backends.PydapDataStore.open(in_url, session=session)
       ds=xa.open_dataset(store)
       print('Succeed accessing MERRA2 OPeNDAP dataset')
    nc_cdate = np.datetime64('%s-%s-%sT%s:30:00'%(yy,mm,dd,hh))
    ds=ds.sel(time=nc_cdate)

    if (area!='Glb'):
       ds=ds.sel(lat=slice(minlat,maxlat),lon=slice(minlon,maxlon))

    outname='%s/%s_%s.%s.png' %(outputpath,area,pltvar,date)
    
    pltdata=ds[pltvar]
    fig,ax=setupax_2dmap(cornerll,area,proj,lbsize=16.)
    set_size(axe_w,axe_h,b=0.13,l=0.05,r=0.95,t=0.95)
    cn=ax.contourf(pltdata.lon,pltdata.lat,pltdata,levels=cnlvs,cmap=clrmap,norm=aer_norm,extend='both')
    ax.set_title(date,loc='left')
    plt.colorbar(cn,ax=ax,orientation='horizontal',ticks=cnlvs[::tkfreq],
                 fraction=0.045,aspect=40,pad=0.08,label=cblb)
    print(outname)
    fig.savefig(outname,dpi=quality)
    plt.close()
 
    dates_count+=1
