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

if (machine=='S4'):
   inputpath='/data/users/swei/common/MERRA2_L64'

sdate=2020060106
edate=2020060106
hint=6
pltall=0 # 0: total only 1: sub species included

varlst=['DU001','DU002','DU003','DU004','DU005',
        'SS001','SS002','SS003','SS004','SS005',
        'OCPHILIC','OCPHOBIC','BCPHILIC','BCPHOBIC',
        'SO4']

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
    c_dtobj=datetime(int(yy),int(mm),int(dd),int(hh))

    date_m3=str(ndate(-3,date))
    m3yy=date_m3[:4] ; m3mm=date_m3[4:6] ; m3dd=date_m3[6:8] ; m3hh=date_m3[8:10]
    m3_dtobj=datetime(int(m3yy),int(m3mm),int(m3dd),int(m3hh))
   
    date_p3=str(ndate( 3,date))
    p3yy=date_p3[:4] ; p3mm=date_p3[4:6] ; p3dd=date_p3[6:8] ; p3hh=date_p3[8:10]
    p3_dtobj=datetime(int(p3yy),int(p3mm),int(p3dd),int(p3hh))

# MERRA2_AER3D_FV3L64.2020060609.nc
    infile_m3=inputpath+'/'+yy+'/'+mm+'/MERRA2_AER3D_FV3L64.'+date_m3+'.nc'
    infile=inputpath+'/'+yy+'/'+mm+'/MERRA2_AER3D_FV3L64.'+date+'.nc'
    infile_p3=inputpath+'/'+yy+'/'+mm+'/MERRA2_AER3D_FV3L64.'+date_p3+'.nc'
   
    ds=xa.open_dataset(infile)
    m3=xa.open_dataset(infile_m3)
    p3=xa.open_dataset(infile_p3)

    an_hour=timedelta(hours=1)
    m3_to_c_dates=pd.date_range(start=m3_dtobj, end=c_dtobj, freq=an_hour)
    c_to_p3_dates=pd.date_range(start=c_dtobj, end=p3_dtobj, freq=an_hour)

    tmpds1=ds.copy()
    tmpds1=tmpds1.assign_coords({'newtime':(['time'],m3_to_c_dates[1].to_datetime64())})
    tmpds1=tmpds1.swap_dims({'time':'newtime'})
    tmpds1=tmpds1.drop_dim('time')

    #tmpds2=m3*1/3+ds*2/3
    #tmpds2.assign(time=m3_to_c_dates[2])
    



 
    dates_count+=1
