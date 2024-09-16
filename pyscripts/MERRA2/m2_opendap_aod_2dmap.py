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
    rootarch='/data/users/swei/ResearchData/Prospectus/AeroObsStats/nc_diag'
    rootpath='/data/users/swei/AlbanyWork/Prospectus/Experiments/HazyDA/Images'
    rootgit='/home/swei/research'
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
tlsize=12 ; txsize=12
mpl.rc('axes', titlesize=tlsize,labelsize=txsize)
mpl.rc('xtick',labelsize=txsize)
mpl.rc('ytick',labelsize=txsize)
mpl.rc('legend',fontsize='small')
fsave=1 ; ffmt='png' ; ptsize=4
axe_w=6 ; axe_h=3 ; quality=300
minussign=u'\u2212'

#inputpath='/glade/work/dfgrogan/UFS/WM_DTAER/AER'
m2_url = 'https://goldsmr4.gesdisc.eosdis.nasa.gov/opendap/MERRA2/M2T1NXAER.5.12.4'
#/2020/09/MERRA2_401.tavg1_2d_aer_Nx.20200901.nc4'
outputpath=rootpath+'/Dataset/MERRA-2/2dmap/AOD'
if ( not os.path.exists(outputpath) ):
    os.makedirs(outputpath)

sdate=2020082200
edate=2020092118
hint=6
pltvar='TOTEXTTAU'
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

date1 = datetime(syy,smm,sdd,shh,30)
date2 = datetime(eyy,emm,edd,ehh,30)
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
fileslist=[]
for date in dlist:
    yy=date[:4] ; mm=date[4:6] ; dd=date[6:8] ; hh=date[8:10]
    pdy=date[:8]
    if (yy=='2020' and mm=='09'):
       m2ind='401'
    else:
       m2ind='400'

#/2020/09/MERRA2_401.tavg1_2d_aer_Nx.20200901.nc4'
    if (pdy!=p_pdy):
       p_pdy=pdy 
       fileslist.append('{}/{}/{}/MERRA2_{}.{}.{}{}{}.nc4'.format(m2_url,yy,mm,m2ind,m2tag, yy, mm, dd))
       
    dates_count+=1

print('Created the list of MERRA-2 files',flush=1)

ds=xa.open_mfdataset(fileslist)
print('Succeed accessing MERRA2 OPeNDAP dataset',flush=1)
ds=ds.sel(time=dates)
if (area!='Glb'):
   ds=ds.sel(lat=slice(minlat,maxlat),lon=slice(minlon,maxlon))

outname='%s/%s_%s.%s_%s.png' %(outputpath,area,pltvar,sdate,edate)

pltdata=ds[pltvar].mean(dim='time')
cblb=ds[pltvar].standard_name
fig,ax,gl=setupax_2dmap(cornerll,area,proj,lbsize=txsize)
set_size(axe_w,axe_h,b=0.13,l=0.05,r=0.95,t=0.95)
cn=ax.contourf(pltdata.lon,pltdata.lat,pltdata,levels=cnlvs,cmap=clrmap,norm=aer_norm,extend='both')
#ax.set_title(loc='left')
plt.colorbar(cn,ax=ax,orientation='horizontal',ticks=cnlvs[::tkfreq],
             fraction=0.045,aspect=40,pad=0.08,label=cblb)
print(outname,flush=1)
fig.savefig(outname,dpi=quality)
plt.close()

