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
import setuparea as setarea
from plot_utils import setupax_2dmap, plt_x2y, set_size
from utils import ndate,setup_cmap
from datetime import datetime, timedelta
from gsi_ncdiag import read_rad_ncdiag

tlsize=12 ; lbsize=10
mpl.rc('axes', titlesize=tlsize,labelsize=lbsize)
mpl.rc('xtick',labelsize=lbsize)
mpl.rc('ytick',labelsize=lbsize)
mpl.rc('legend',fontsize='medium')
fsave=1 ; ffmt='png' ; ptsize=4
axe_w=3 ; axe_h=3 ; quality=300
minussign=u'\u2212'

# Plotting setup
sdate=2020061000
edate=2020092118
hint=6
exp='observer_v2qc'
saved_suffix='v2qc'
sensor='iasi_metop-a'
spectral_range=slice(700,1300)
loop='ges' #ges,anl
usebc=0
pltline=1 # plot 2d histogram
#chkwvn_list=[906.25]
chkwvn_list=[906.25,943.25,962.5,1096.0,]

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
savedir=outpath+'/DiagFiles/rad/observer/aelines_'+saved_suffix
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
        print('Processing Radfile: %s' %(raddfile),flush=1)
        dates_count+=1
    else:
        print('%s is not existing'%(raddfile),flush=1)
        continue
    
    tmpds=read_rad_ncdiag(infile,chkwvn=chkwvn_list,cal_ae=1,get_water_frac=1)
    
    if (date==str(sdate)):
        ds_all=tmpds
    else:
        ds_all=xa.concat((ds_all,tmpds),dim='obsloc')

print('Obs. counts: %i' %(ds_all.obsloc.size), flush=1)

#for chkwvn in [962.5]:
for chkwvn in chkwvn_list:
    ds_chk=ds_all.sel(wavenumber=chkwvn)
    tb_sim=ds_chk.tb_sim
    tb_clr=ds_chk.tb_clr
    tb_obs=ds_chk.tb_obs
    if (usebc):
       omb=ds_chk.omb_bc
       omb_clr=(tb_obs-tb_clr)-(omb-tb_obs+tb_sim)
    else:
       omb=ds_chk.omb_nbc
       omb_clr=tb_obs-tb_clr

    aereff_fg=tb_sim-tb_clr
    aereff_obs=tb_obs-tb_clr
    aereff=0.5*abs(aereff_fg)+0.5*abs(aereff_obs)

    water_msk=(ds_chk.water_frac>0.999)
    good_msk=(ds_chk.qcflag==0.)
    gross_msk=(ds_chk.qcflag==3.)
    cld_msk=(ds_chk.qcflag==7.)
    tzr_msk=(ds_chk.qcflag==10.)
    sfcir_msk=(ds_chk.qcflag==53.)
    bust_msk=(ds_chk.qcflag==55.)
    aercld_msk=(ds_chk.qcflag==57.)
    if (saved_suffix=='v1qc'):
       bustqc=(abs(omb)>3.)&(abs(omb)>1.8*ds_chk.Ae)
       tmpmsk=(aercld_msk)&(~bustqc)
       omb=xa.where(tmpmsk,omb_clr,omb)
       aer_msk=(ds_chk.qcflag==13.)|(tmpmsk)
    else:
       aer_msk=(ds_chk.qcflag==13.)
    passed_msk=(good_msk)|(aer_msk)
    # ori_msk=(good_msk)|(aer_msk)|(bust_msk)|(tzr_msk)

    pltmask=passed_msk
    
    qced_flg=np.zeros_like(passed_msk)
    qced_flg[passed_msk]=1.
    wrkds=xa.Dataset({'ae_fg':(['obsloc'],abs(aereff_fg).data),
                      'ae_obs':(['obsloc'],abs(aereff_obs).data),
                      'ae_sym':(['obsloc'],aereff.data),
                      'pltflg':(['obsloc'],qced_flg),
                      'omb':(['obsloc'],omb.data),
                      },coords={'obsloc':np.arange(aereff.size)})

    if (pltline):
        counts=np.count_nonzero(pltmask)
        print('%.2f cm-1: %i' %(chkwvn,counts))

        x_label='Aerosol effect [K]'
        if (usebc):
           y_label='Mean OMFs w/ BC [K]'
        else:
           y_label='Mean OMFs w/o BC [K]'

        binsize=0.1
        halfbin=0.5*binsize
        hist_x_edge=np.arange(-1*halfbin,10.+binsize,binsize)
        bin_center=(hist_x_edge+halfbin)[:-1]

        tistr=('%s (%.2f $cm^{-1}$)' %(sensor,chkwvn))

        df0=wrkds.to_dataframe()
        df0=df0.reset_index()
        filter=((df0['pltflg']==1.))
        wrkdf=df0.loc[filter,:]
        
        wrkdf['ae_fg_bin']=pd.cut(wrkdf['ae_fg'],bins=hist_x_edge,labels=bin_center)
        wrkdf['ae_obs_bin']=pd.cut(wrkdf['ae_obs'],bins=hist_x_edge,labels=bin_center)
        wrkdf['ae_sym_bin']=pd.cut(wrkdf['ae_sym'],bins=hist_x_edge,labels=bin_center)
      
        b=0
        for bin in ['ae_fg_bin','ae_obs_bin','ae_sym_bin']:
            tmpgrp=wrkdf.groupby([bin]).agg({'omb':['mean']}).reset_index()
            tmpgrp=tmpgrp.droplevel(1,axis=1)
            newcol=bin[:-4]+'_omb'
            tmpgrp=tmpgrp.rename(columns={'omb':newcol})
            if b==0:
               pltgrp=tmpgrp
            else:
               pltgrp[newcol]=tmpgrp[newcol]
            b+=1

        fig,ax=plt.subplots()
        set_size(axe_w,axe_h,l=0.15,b=0.15)
        ax.set_prop_cycle(color=['r','b','k'])
        pltgrp.plot(x='ae_fg_bin',ax=ax,alpha=0.8,ms=2,zorder=4)
        #pltgrp[['ae_fg_omb','ae_obs_omb','ae_sym_omb']].plot(x='ae_fg_bin',ax=ax,alpha=0.8,ms=2,zorder=4)
        ax.legend(['Ae_FG','Ae_OBS','Ae_Sym'])
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.set_title(tistr)
        xmin,xmax=ax.get_xbound()
        #ax.hlines(0.,xmin,xmax,colors='k',lw=0.8,ls='--',zorder=3)
        ax.set_xlim(xmin,xmax)

        if (fsave):
            fname=('%s/Lines_%s_%s_%.2f_%s.%s_%s.%s'
                    %(savedir,area,sensor,chkwvn,bcflg,sdate,edate,ffmt))
            print(fname,flush=1)
            fig.savefig(fname,dpi=quality)
            plt.close()
