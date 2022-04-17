# -*- coding: utf-8 -*-
"""
Created on Tue Jun 29 16:06:23 2021

@author: ck102

Aerosol detection based on CADS 3.1 from NWP SAF

"""
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
    rootpath='/data/users/swei'
    rootarch='/scratch/users/swei/ncdiag'
    rootgit='/home/swei/research'
elif (machine=='Hera'):
    rootpath='/scratch2/BMC/gsd-fv3-dev/Shih-wei.Wei'
    rootarch='/scratch2/BMC/gsd-fv3-dev/Shih-wei.Wei/ResearchData'
    rootgit='/home/Shih-wei.Wei/research'
elif (machine=='Cheyenne'):
    rootpath='/glade/work/swei/output/images'
    rootarch='/scratch2/BMC/gsd-fv3-dev/Shih-wei.Wei/ResearchData'
    rootgit='/glade/u/home/swei/research'
sys.path.append(rootgit+'/pyscripts/functions')
import numpy as np
import xarray as xa
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mpcrs
import setuparea as setarea
from plot_utils import setupax_2dmap, plt_x2y, set_size
from utils import ndate,setup_cmap
from datetime import datetime, timedelta

tlsize=12 ; lbsize=10
mpl.rc('axes', titlesize=tlsize,labelsize=lbsize)
mpl.rc('xtick',labelsize=lbsize)
mpl.rc('ytick',labelsize=lbsize)
mpl.rc('legend',fontsize='large')
fsave=1 ; ffmt='png' ; ptsize=4
axe_w=3 ; axe_h=3 ; quality=300
minussign=u'\u2212'

# Plotting setup
sdate=2020061000
edate=2020071018
hint=6
exp='hazyda_aero'
sensor='iasi_metop-a'
spectral_range=slice(700,1300)
loop='ges' #ges,anl
usebc=0
plthist2d=1 # plot 2d histogram

area='Glb'
minlon, maxlon, minlat, maxlat, crosszero, cyclic=setarea.setarea(area)
print(minlat,maxlat,minlon,maxlon,crosszero,cyclic)
if (area=='Glb'):
   minlon=-180. ; maxlon=180.
cornll=[minlat,maxlat,minlon,maxlon]

cbori='vertical' #vertical, horizontal
if (cbori=='vertical'):
   cb_frac=0.025
   cb_pad=0.06
elif (cbori=='horizontal'):
   cb_frac=0.04
   cb_pad=0.1

if (usebc):
   bcflg='bc'
else:
   bcflg='nobc'

# Data path setup
outpath=rootpath+'/AlbanyWork/Prospectus/Experiments/HazyDA/Images'
archdir=rootarch+'/'+exp
savedir=outpath+'/DiagFiles/rad/hist2d'
if ( not os.path.exists(savedir) ):
    os.makedirs(savedir)

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
    raddfile='diag_'+sensor+'_'+loop+'.'+str(date)+'.nc4'
    infile=archdir+'/'+str(date)+'/'+raddfile

    if (os.path.exists(infile)):
        print('Processing Radfile: %s' %(raddfile))
        ds=xa.open_dataset(infile)
        npts=int(ds.nobs.size/ds.nchans.size)
        nchs=ds.nchans.size
        ds=ds.assign_coords(nuchan=('wavenumber',ds.wavenumber.data))
        ds=ds.swap_dims({"nchans":"wavenumber"}) #replace the dimension of channel by channel indices
        wavelength=1e+04/ds.wavenumber
        chkwvn_list=ds.wavenumber.sel(wavenumber=spectral_range)[ds.use_flag.sel(wavenumber=spectral_range)==1]
        dates_count+=1
    else:
        print('%s is not existing'%(raddfile))
        continue
    
    # Observation lat/lon
    rlat=np.reshape(ds.Latitude.values,(npts,nchs))
    rlon=np.reshape(ds.Longitude.values,(npts,nchs))
    qcflags=np.reshape(ds.QC_Flag.values,(npts,nchs))
    #obs=np.reshape(ds.Observation.values,(npts,nchs))
    sim=np.reshape(ds.Simulated_Tb.values,(npts,nchs))
    clr=np.reshape(ds.Clearsky_Tb.values,(npts,nchs))
    omb_bc=np.reshape(ds.Obs_Minus_Forecast_adjusted.values,(npts,nchs))
    omb_nbc=np.reshape(ds.Obs_Minus_Forecast_unadjusted.values,(npts,nchs))
    obs=omb_nbc+sim
    tmpds=xa.Dataset({'rlon':(['obsloc'],rlon[:,0]),
                      'rlat':(['obsloc'],rlat[:,0]),
                      'qcflag':(['obsloc','wavenumber'],qcflags),
                      'tb_obs':(['obsloc','wavenumber'],obs),
                      'tb_sim':(['obsloc','wavenumber'],sim),
                      'tb_clr':(['obsloc','wavenumber'],clr),
                      'omb_bc':(['obsloc','wavenumber'],omb_bc),
                      'omb_nbc':(['obsloc','wavenumber'],omb_nbc),
                      },
                     coords={'obsloc':np.arange(npts),
                             'wavenumber':ds.wavenumber.values})
    
    tmpds=tmpds.sel(wavenumber=chkwvn_list)
    
    if (date==str(sdate)):
        ds_all=tmpds
    else:
        ds_all=xa.concat((ds_all,tmpds),dim='obsloc')

#for chkwvn in [962.5]:
for chkwvn in chkwvn_list:
    ds_chk=ds_all.sel(wavenumber=chkwvn)
    tb_sim=ds_chk.tb_sim
    tb_clr=ds_chk.tb_clr
    tb_obs=ds_chk.tb_obs
    if (usebc):
       omb=ds_chk.omb_bc
    else:
       omb=ds_chk.omb_nbc

    aereff_fg=tb_sim-tb_clr
    aereff_obs=tb_obs-tb_clr
    aereff=0.5*abs(aereff_fg)+0.5*abs(aereff_obs)
    
    qc0_msk=(ds_chk.qcflag==0.)
    qc7_msk=(ds_chk.qcflag==7.)
    qc13_msk=(ds_chk.qcflag==13.)

    pltmask=(qc0_msk)|(qc13_msk)
    if (plthist2d):
        for pltset in [2]:
            if (pltset==1):
                pltda_x=tb_sim[pltmask==1]
                pltda_y=tb_obs[pltmask==1]
                x_label='Background BT [K]'
                y_label='Observation BT [K]'
                hist_x_edge=np.arange(200,300.5,0.5)
                hist_y_edge=hist_x_edge
            elif (pltset==2):
                pltda_x=aereff[pltmask==1]
                pltda_y=omb[pltmask==1]
                x_label='Aerosol effect [K]'
                if (usebc):
                   y_label='OMB w/ BC [K]'
                else:
                   y_label='OMB w/o BC [K]'
                hist_x_edge=np.arange(0.,40.5,0.5)
                hist_y_edge=np.arange(-29.75,30.25,0.5)
         
            counts=np.count_nonzero(pltmask)
        # weights_arr=np.zeros((hist_x_edge.size),dtype='float')
        # weights_arr[:]=np.count_nonzero(pltmask)
        
        # print('%s used %i of %i' %(sensor,usedchidx.size,ds1.channel.size))
        
        #tistr=('%s %s %s %s: AOD>%.2f %s>%.2f (%.2fµm)'
        #       %(imidx,tlstr,sensor,area,aodmin,aersp,aerfracmin,wavelength[chkwvlidx]))
            tistr=('%s (%.2f $cm^{-1}$)' %(sensor,chkwvn))
            
            cbtcks=np.arange(0,5.1,0.1)
            lvs=np.power(10,cbtcks)
            clridx=[0]
            for idx in np.linspace(2,128,cbtcks.size-1):
                clridx.append(int(idx))
            clrmap=setup_cmap('MPL_jet',clridx)
            norm = mpcrs.BoundaryNorm(lvs,len(clridx)+1,extend='both')
            
            fig=plt.figure()
            ax=plt.subplot()
            set_size(axe_w,axe_h,l=0.15)
            ax.set_title(tistr,loc='left')
            # hdata,xedgs,yedgs,im=ax.hist2d(pltda_x,pltda_y,[hist_x_edge,hist_y_edge],cmap=clrmap,norm=norm)
            hdata,xedgs,yedgs,im=ax.hist2d(pltda_x,pltda_y,
                                           [hist_x_edge,hist_y_edge],
                                           density=0,cmap=clrmap,
                                           norm=norm)
            # cbar=plt.colorbar(im,orientation=cbori,fraction=0.05,ticks=lvs,label=r'Probability[%]',pad=0.02)
            cbar=plt.colorbar(im,orientation=cbori,ticks=lvs[::10],fraction=cb_frac,pad=cb_pad,aspect=30)
            if (cbori=='horizontal'):
                cbar.ax.set_xticklabels(cbtcks[::10])
            else:
                cbar.ax.set_yticklabels(cbtcks[::10])
                
            #if (pltset==2):
            #    bndry_line_ny=np.arange(-30,-2)
            #    bndry_line_nx=bndry_line_ny/1.8
            #    bndry_line_py=np.arange(3,31)
            #    bndry_line_px=bndry_line_py/1.8
            #    plt.plot(abs(bndry_line_nx),bndry_line_ny,'k',linestyle='dashed')
            #    plt.plot(bndry_line_px,bndry_line_py,'k',linestyle='dashed')
                
            ax.set_xlabel(x_label)
            ax.set_ylabel(y_label)
            
            if (fsave):
                fname=('%s/2DPDF_%s_%s_%.2f_%s.set%s.%s_%s.%s'
                        %(savedir,area,sensor,chkwvn,bcflg,pltset,sdate,edate,ffmt))
                print(fname)
                fig.savefig(fname,dpi=quality)
                plt.close()
