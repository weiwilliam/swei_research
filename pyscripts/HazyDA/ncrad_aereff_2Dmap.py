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
    rootpath='/data/users/swei/AlbanyWork/Prospectus/Experiments/HazyDA/Images'
    rootarch='/data/users/swei/ResearchData/Prospectus/AeroObsStats/nc_diag'
    rootgit='/home/swei/research'
elif (machine=='Hera'):
    rootpath='/scratch2/BMC/gsd-fv3-dev/Shih-wei.Wei'
    rootarch='/scratch2/BMC/gsd-fv3-dev/Shih-wei.Wei/ResearchData'
    rootgit='/home/Shih-wei.Wei/research'
elif (machine=='Cheyenne'):
    rootpath='/glade/work/swei/output/images'
    rootarch='/scratch2/BMC/gsd-fv3-dev/Shih-wei.Wei/ResearchData'
    rootgit='/glade/u/home/swei/research'
import numpy as np
import xarray as xa
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mpcrs
import cartopy.crs as ccrs
from datetime import datetime, timedelta
import setuparea as setarea
from plot_utils import setupax_2dmap, plt_x2y, set_size
from utils import ndate,setup_cmap
from gsi_ncdiag import read_rad_ncdiag

tlsize=12 ; lbsize=10
mpl.rc('axes', titlesize=tlsize,labelsize=lbsize)
mpl.rc('xtick',labelsize=lbsize)
mpl.rc('ytick',labelsize=lbsize)
mpl.rc('legend',fontsize='large')
fsave=1 ; ffmt='png' ; ptsize=2
axe_w=6 ; axe_h=3 ; quality=300

# Projection setting
proj=ccrs.PlateCarree(globe=None)

# Plotting setup
exp='observer_v2qc'
sdate=2020062212
edate=2020062212
aertype='Dust'
hint=6
sensor='iasi_metop-a'
spectral_range=slice(700,1300)
chkwvn=962.50
loop='ges' #ges,anl
#if loop=='anl':
#    tlstr='OMA'
#elif loop=='ges':
#    tlstr='OMF'
plt2d_ae=1  # plot single cycle 2d map ae
plt2d_omb=0  # plot single cycle 2d map omb

area='Glb'
minlon, maxlon, minlat, maxlat, crosszero, cyclic=setarea.setarea(area)
print(minlat,maxlat,minlon,maxlon,crosszero,cyclic)
if (area=='Glb'):
   minlon=-180. ; maxlon=180.
cornll=[minlat,maxlat,minlon,maxlon]

cbori='horizontal' #vertical, horizontal
if (cbori=='vertical'):
   cb_frac=0.025
   cb_pad=0.04
elif (cbori=='horizontal'):
   cb_frac=0.04
   cb_pad=0.1

# Data path setup
archdir=rootarch+'/'+exp
outpath=rootpath+'/DiagFiles/rad/observer'

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

naer=14
dates_count=0
for date in dlist:
    raddfile='diag_'+sensor+'_'+loop+'.'+str(date)+'.nc4'
    infile=archdir+'/'+str(date)+'/'+raddfile

    if (os.path.exists(infile)):
        print('Processing Radfile: %s' %(raddfile))
        dates_count+=1
    else:
        print('%s is not existing'%(raddfile))
        continue
    
    tmpds=read_rad_ncdiag(infile,chkwvn=chkwvn)
    
    tb_sim=tmpds.tb_sim
    tb_clr=tmpds.tb_clr
    tb_obs=tmpds.tb_obs

    omb=tb_obs-tb_sim
    aereff_fg=tb_sim-tb_clr
    aereff_obs=tb_obs-tb_clr
    aereff=0.5*abs(aereff_fg)+0.5*abs(aereff_obs)
   
    all_msk=(~np.isnan(tmpds.obsloc))
    good_msk=(tmpds.qcflag==0.)
    hazy_msk=(tb_sim!=tb_clr)
    #cld_msk=(tmpds.qcflag==7.)
    aer_msk=(tmpds.qcflag==13.)
    aergrs_msk=(tmpds.qcflag==55.)
    aercld_msk=(tmpds.qcflag==57.)

    #pltmask=all_msk
    #pltmask=(good_msk)|(aer_msk)|(aergrs_msk)|(aercld_msk)             
    pltmask=(good_msk)|(hazy_msk)
    
    if (plt2d_ae):
        savedir=outpath+'/2d_Ae'
        if ( not os.path.exists(savedir) ):
            os.makedirs(savedir)
        
        cbtcks=np.arange(0,11,0.5)
        lvs=cbtcks
        clridx=[0]
        for idx in np.linspace(64,128,cbtcks.size-1):
            clridx.append(int(idx))
        clrmap=setup_cmap('MPL_jet',clridx)
        norm = mpcrs.BoundaryNorm(lvs,len(clridx)+1,extend='both')
        
        pltda=aereff[pltmask==1]
        cblabel='Aerosol effect [K]'
       
        tistr=('%s (%.2f $cm^{-1}$)' %(sensor,chkwvn))
        
        fig,ax,gl=setupax_2dmap(cornll,area,proj,12)
        set_size(axe_w,axe_h)
        ax.set_title(tistr,loc='left')
        sc=ax.scatter(tmpds.rlon[pltmask],tmpds.rlat[pltmask],c=pltda,
                   s=ptsize,cmap=clrmap,norm=norm)
        plt.colorbar(sc,orientation=cbori,fraction=cb_frac,
                     pad=cb_pad,ticks=lvs[::2],aspect=40,label=cblabel)
        
        if (fsave):
            fname=('%s/%s_%s_%.2f.%s.%s' %(savedir,area,sensor,chkwvn,str(date),ffmt))
            print(fname)
            fig.savefig(fname,dpi=quality)
            plt.close()

        cbtcks=np.arange(-10,11,1)
        lvs=cbtcks
        clridx=[]
        for idx in np.linspace(2,254,cbtcks.size):
            clridx.append(int(idx))
        clrmap=setup_cmap('BlueYellowRed',clridx)
        norm = mpcrs.BoundaryNorm(lvs,len(clridx)+1,extend='both')

        pltda=aereff_fg[pltmask==1]
        cblabel='Aerosol effect FG [K]'

        tistr=('%s (%.2f $cm^{-1}$)' %(sensor,chkwvn))

        fig,ax,gl=setupax_2dmap(cornll,area,proj,12)
        set_size(axe_w,axe_h)
        ax.set_title(tistr,loc='left')
        sc=ax.scatter(tmpds.rlon[pltmask],tmpds.rlat[pltmask],c=pltda,
                   s=ptsize,cmap=clrmap,norm=norm)
        plt.colorbar(sc,orientation=cbori,fraction=cb_frac,
                     pad=cb_pad,ticks=lvs[::2],aspect=40,label=cblabel)

        if (fsave):
            fname=('%s/%s_%s_AeFG_%.2f.%s.%s' %(savedir,area,sensor,chkwvn,str(date),ffmt))
            print(fname)
            fig.savefig(fname,dpi=quality)
            plt.close()

        pltda=aereff_obs[pltmask==1]
        cblabel='Aerosol effect OBS [K]'

        tistr=('%s (%.2f $cm^{-1}$)' %(sensor,chkwvn))

        fig,ax,gl=setupax_2dmap(cornll,area,proj,12)
        set_size(axe_w,axe_h)
        ax.set_title(tistr,loc='left')
        sc=ax.scatter(tmpds.rlon[pltmask],tmpds.rlat[pltmask],c=pltda,
                   s=ptsize,cmap=clrmap,norm=norm)
        plt.colorbar(sc,orientation=cbori,fraction=cb_frac,
                     pad=cb_pad,ticks=lvs[::2],aspect=40,label=cblabel)

        if (fsave):
            fname=('%s/%s_%s_AeOBS_%.2f.%s.%s' %(savedir,area,sensor,chkwvn,str(date),ffmt))
            print(fname)
            fig.savefig(fname,dpi=quality)
            plt.close()

    if (plt2d_omb):
        savedir=outpath+'/2d_OMB'
        if ( not os.path.exists(savedir) ):
            os.makedirs(savedir)
        
        cbtcks=np.arange(-10,11,1)
        lvs=cbtcks
        clridx=[]
        for idx in np.linspace(2,254,cbtcks.size):
            clridx.append(int(idx))
        clrmap=setup_cmap('BlueYellowRed',clridx)
        norm = mpcrs.BoundaryNorm(lvs,len(clridx)+1,extend='both')
        
        pltda=omb[pltmask==1]
        cblabel='Aerosol effect [K]'
       
        tistr=('%s (%.2f $cm^{-1}$)' %(sensor,chkwvn))
        
        fig,ax,gl=setupax_2dmap(cornll,area,proj,12)
        set_size(axe_w,axe_h)
        ax.set_title(tistr,loc='left')
        sc=ax.scatter(tmpds.rlon[pltmask],tmpds.rlat[pltmask],c=pltda,
                   s=ptsize,cmap=clrmap,norm=norm)
        plt.colorbar(sc,orientation=cbori,fraction=cb_frac,pad=cb_pad,ticks=lvs,aspect=40)
        
        if (fsave):
            fname=('%s/%s_%s_%.2f.%s.%s' %(savedir,area,sensor,chkwvn,str(date),ffmt))
            print(fname)
            fig.savefig(fname,dpi=quality)
            plt.close()
