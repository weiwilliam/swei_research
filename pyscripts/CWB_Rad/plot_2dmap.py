#!/usr/bin/env python3
import sys, os, platform
import numpy as np
import pandas as pd
import xarray as xa
from datetime import datetime, timedelta
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mpcrs
import cartopy.crs as ccrs
import cartopy.feature as cft

from utils import setup_cmap, ndate
from plot_utils import setupax_2dmap,set_size
from TCo_latlon_def import define_latlon
import setuparea as setarea

# Plotting setup
txsize=12
mpl.rc('axes',titlesize=12,labelsize=12)
mpl.rc('xtick',labelsize=12)
mpl.rc('ytick',labelsize=12)
mpl.rc('legend',fontsize='large')
axe_w=6; axe_h=4
quality=300

inputpath='/glade/work/dfgrogan/UFS/WM_DTAER/AER'

sdate=2020082200
edate=2020092118
hint=6
pltvar='AODANA'
area='NAmer'
pltave=0 # 0: single cycle only 1: time average
tkfreq=1

#outputpath=rootpath+'/Dataset/MERRA-2/2dmap/AOD_new/'+area
#if ( not os.path.exists(outputpath) ):
#    os.makedirs(outputpath)
#
#cnlvs=np.array((0., 0.05, 0.1, 0.15, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 2.5))
#clridx=np.array((0,2,4,7,8,9,10,12,14,15,16,17,18))
#clrmap=setup_cmap('precip3_16lev',clridx)
#aer_norm = mpcrs.BoundaryNorm(cnlvs,len(clridx),extend='both')
#cblb='Analysis AOD'

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

date1=pd.to_datetime(sdate,format="%Y%m%d%H")
date2=pd.to_datetime(edate,format="%Y%m%d%H")
delta = timedelta(hours=hint)
dates = pd.date_range(start=date1, end=date2, freq=delta)

lat,lon=define_latlon(jcap=383)

sys.exit(0)

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
    infile=inputpath+'/MERRA2_'+m2ind+'.'+m2tag+'.'+pdy+'_'+hh+'Z.nc4'
    ds=xa.open_dataset(infile)
    tmpvar=ds[pltvar]
    if (area!='Glb'):
       tmpvar=tmpvar.sel(lon=slice(minlon,maxlon),lat=slice(minlat,maxlat))

    outname='%s/%s_%s.%s.png' %(outputpath,area,pltvar,date)
    
    pltdata=tmpvar.sel(time=ds.time[0])
    fig,ax,gl=setupax_2dmap(cornerll,area,proj,lbsize=txsize)
    set_size(axe_w,axe_h,b=0.13,l=0.05,r=0.95,t=0.95)
    cn=ax.contourf(pltdata.lon,pltdata.lat,pltdata,levels=cnlvs,cmap=clrmap,norm=aer_norm,extend='both')
    ax.set_title(date,loc='left')
    plt.colorbar(cn,ax=ax,orientation='horizontal',ticks=cnlvs[::tkfreq],
                 fraction=0.045,aspect=40,pad=0.08,label=cblb)
    if (area=='NAmer'):
       ax.add_feature(cft.STATES,zorder=2)
    print(outname)
    fig.savefig(outname,dpi=quality)
    plt.close()

    if (pltave):
       if (dates_count==0):
          var=tmpvar
       else:
          var=xa.concat((var,tmpvar),dim='time')

       if (date==str(edate)):
          outname='%s/%s_%s.mean.%s_%s.png' %(outputpath,area,pltvar,sdate,edate)
          title_str='%s-%s'%(sdate,edate)

          pltdata=var.mean(dim='time')
          fig,ax,gl=setupax_2dmap(cornerll,area,proj,lbsize=txsize)
          set_size(axe_w,axe_h,b=0.13,l=0.05,r=0.95,t=0.95)
          cn=ax.contourf(pltdata.lon,pltdata.lat,pltdata,levels=cnlvs,cmap=clrmap,norm=aer_norm,extend='both')
          ax.set_title(title_str,loc='left')
          plt.colorbar(cn,ax=ax,orientation='horizontal',ticks=cnlvs[::tkfreq],
                       fraction=0.045,aspect=40,pad=0.08,label=cblb)
          if (area=='NAmer'):
             ax.add_feature(cft.STATES,zorder=2)
          print(outname)
          fig.savefig(outname,dpi=quality)
          plt.close()
 
    dates_count+=1
