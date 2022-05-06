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
import matplotlib.pyplot as plt

# Plotting setup
expname='hazyda_aero'
cdate=2020062212
hint=6
sensor='iasi_metop-a'
selwvn=962.5
loop='ges' #ges,anl
degres=1

raddfile='diag_'+sensor+'_'+loop+'.'+str(cdate)+'.nc4'
#infile1=inputpath+'/'+expname+'/'+str(cdate)+'/'+raddfile
tmppath='/data/users/swei/FTPdir/'
# tmppath='/data/users/swei/Experiments/testing/OUTPUT/bc_check/2020062212'
# tmppath='F:/ResearchData/Prospectus/HazyDA'
#tmppath='/scratch/users/swei/ncdiag/hazyda_aero/'+str(cdate)
infile1=tmppath+'/'+raddfile

if (os.path.exists(infile1)):
    print('Processing Radfile: %s' %(raddfile))
    ds1=xa.open_dataset(infile1)
    npts=int(ds1.nobs.size/ds1.nchans.size)
    nchs=ds1.nchans.size
    chkwvn_list=ds1.wavenumber[ds1.use_flag==1]
    #ds1=ds1.swap_dims({"nchans":"wavenumber"}) #replace the dimension of channel by channel indices
    aerdiag=0
    #if ('aero_load' in ds1.variables):
        #naer=ds1.aero_frac_arr_dim.size
        #aerdiag=1
    
    rlat1=np.reshape(ds1.Latitude.values,(npts,nchs))
    rlon1=np.reshape(ds1.Longitude.values,(npts,nchs))
    rlon1=(rlon1+180)%360-180
    qcflags=np.reshape(ds1.QC_Flag.values,(npts,nchs))
    sim1=np.reshape(ds1.Simulated_Tb.values,(npts,nchs))
    clr1=np.reshape(ds1.Clearsky_Tb.values,(npts,nchs))
    sim_nbc1=np.reshape(ds1.Obs_Minus_Forecast_unadjusted.values,(npts,nchs))
    obs1=sim_nbc1+sim1
    bcemiss=np.reshape(ds1.BC_Emissivity.values,(npts,nchs))
    #bcpredemiss=np.reshape(ds1.BCPred_Emissivity.values,(npts,nchs))
    if (aerdiag):
        aero_load=np.reshape(ds1.aero_load.values,(npts,nchs))#[:,0]
        aero_frac=np.reshape(ds1.aero_frac.values,(npts,nchs,naer))#[:,0,:]
        aero_name=ds1.aero_name.split()

                      #'bcpred_emiss':(['obsloc','wavenumber'],bcpredemiss),
    tmpds=xa.Dataset({'rlon':(['obsloc'],rlon1[:,0]),
                      'rlat':(['obsloc'],rlat1[:,0]),
                      'qcflag':(['obsloc','wavenumber'],qcflags),
                      'tb_obs':(['obsloc','wavenumber'],obs1),
                      'tb_sim':(['obsloc','wavenumber'],sim1),
                      'tb_clr':(['obsloc','wavenumber'],clr1),
                      'bc_emiss':(['obsloc','wavenumber'],bcemiss),
                      },
                     coords={'obsloc':np.arange(npts),
                             'wavenumber':ds1.wavenumber.values})

    if (aerdiag):
       tmpds=tmpds.assign({'aero_load':(['obsloc','wavenumber'],aero_load)})
       tmpds=tmpds.assign_coords({'naer':aero_name})
       tmpds=tmpds.assign({'aero_frac':(['obsloc','wavenumber','naer'],aero_frac)})

    #tmpds=tmpds.sel(wavenumber=chkwvn_list.data) 
else:
    print('No such file: %s' %(infile1))

