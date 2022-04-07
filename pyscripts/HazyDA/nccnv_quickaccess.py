# -*- coding: utf-8 -*-
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
        rootarch='/scratch/users/swei/ncdiag'
        rootgit='/home/swei/research'
        machine='S4'
sys.path.append(rootgit+'/pyscripts/functions')
import numpy as np
import xarray as xa
import pandas as pd
from plot_utils import setupax_2dmap, plt_x2y, set_size
from utils import ndate,setup_cmap
from datetime import datetime, timedelta
import setuparea as setarea
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mpcrs
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

# Plotting setup
cdate=2018071418
hint=6
cnvvar='t'
loop='ges' #ges,anl
area='NYS'
minlon, maxlon, minlat, maxlat, crosszero, cyclic=setarea.setarea(area)
print(area,minlat,maxlat,minlon,maxlon,crosszero)
if (area=='Glb'):
   minlon=-180. ; maxlon=180.
else:
   minlon=(minlon+180)%360-180
   maxlon=(maxlon+180)%360-180
cornerll=[minlat,maxlat,minlon,maxlon]

tmppath='/data/users/swei/FTPdir/'
tmppath='C:/Users/ck102/Downloads'
#tmppath='/data/users/swei/archive/hazyda_aero/2020060112'
cnvdfile='diag_conv_'+cnvvar+'_'+loop+'.'+str(cdate)+'.nc4'
#infile1=inputpath+'/'+expname+'/'+str(cdate)+'/'+raddfile
infile1=tmppath+'/'+cnvdfile

if (os.path.exists(infile1)):
    print('Processing Cnvfile: %s' %(cnvdfile))
    ds1=xa.open_dataset(infile1)
    
    npts=ds1.nobs.size
    rlat=ds1.Latitude.data
    rlon=ds1.Longitude.data
    rlon=(rlon+180)%360-180
    sta_id=ds1.Station_ID.data
    obstype=ds1.Observation_Type.data
    sta_elev=ds1.Station_Elevation.data
    qcflags=ds1.Analysis_Use_Flag.data
    errinv=ds1.Errinv_Final.data
    obs=ds1.Observation.data
    omb_bc=ds1.Obs_Minus_Forecast_adjusted.data
    omb_nbc=ds1.Obs_Minus_Forecast_unadjusted.data

    tmpds=xa.Dataset({'rlon':(['obsloc'],rlon),
                      'rlat':(['obsloc'],rlat),
                      'qcflag':(['obsloc'],qcflags),
                      'obs':(['obsloc'],obs),
                      'omb_bc':(['obsloc'],omb_bc),
                      'omb_nbc':(['obsloc'],omb_nbc),
                      'sta_id':(['obsloc'],sta_id),
                      'obstype':(['obsloc'],obstype),
                      },
                     coords={'obsloc':np.arange(npts)})
     
else:
    print('No such file: %s' %(infile1))

pltmsk=tmpds.obstype==188
fig,ax,gl=setupax_2dmap(cornerll,area,proj,lbsize=txsize)
set_size(axe_w,axe_h,b=0.15,l=0.05,r=0.95)
ax.scatter(tmpds.rlon[pltmsk],tmpds.rlat[pltmsk],c='r')

