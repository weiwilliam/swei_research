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

# Plotting setup
cdate=2022030618
hint=6
cnvvar='t'
loop='ges' #ges,anl

tmppath='/data/users/swei/FTPdir/'
#tmppath='/data/users/swei/archive/hazyda_aero/2020060112'
cnvdfile='diag_conv_'+cnvvar+'_'+loop+'.'+str(cdate)+'.nc4'
#infile1=inputpath+'/'+expname+'/'+str(cdate)+'/'+raddfile
infile1=tmppath+'/'+cnvdfile

if (os.path.exists(infile1)):
    print('Processing Cnvfile: %s' %(cnvdfile))
    ds1=xa.open_dataset(infile1)
    
    npts=ds1.nobs.size
    rlat=ds1.Latitude.values
    rlon=ds1.Longitude.values
    rlon=(rlon+180)%360-180
    sta_id=ds1.Station_ID.values
    obstype=ds1.Observation_Type.values
    sta_elev=ds1.Station_Elevation.values
    qcflags=ds1.Analysis_Use_Flag.values
    errinv=ds1.Errinv_Final.values
    obs=ds1.Observation.values
    omb_bc=ds1.Obs_Minus_Forecast_adjusted
    omb_nbc=ds1.Obs_Minus_Forecast_unadjusted

    tmpds=xa.Dataset({'rlon':(['obsloc'],rlon),
                      'rlat':(['obsloc'],rlat),
                      'qcflag':(['obsloc'],qcflags),
                      'obs':(['obsloc'],obs),
                      'omb_bc':(['obsloc'],sim1),
                      'omb_nbc':(['obsloc'],clr1)
                      },
                     coords={'obsloc':np.arange(npts)})
     
else:
    print('No such file: %s' %(infile1))
