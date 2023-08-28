#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 18 12:41:00 2018

@author: weiwilliam
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
    rootarch='/data/users/swei/archive/nc_DiagFiles'
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
import matplotlib.dates as mdates
from matplotlib.dates import (DAILY, DateFormatter,
                              rrulewrapper, RRuleLocator)
import setuparea as setarea
from plot_utils import setupax_2dmap, set_size
from utils import ndate,setup_cmap
from datetime import datetime, timedelta
from opendap_m2 import opendap_m2_aod
import scipy.stats

tlsize=12 ; lbsize=10
mpl.rc('axes', titlesize=tlsize,labelsize=lbsize)
mpl.rc('xtick',labelsize=lbsize)
mpl.rc('ytick',labelsize=lbsize)
mpl.rc('legend',fontsize='small')
fsave=1 ; ffmt='png' ; ptsize=4
axe_w=3 ; axe_h=3 ; quality=300
diagsuffix='nc4'
minussign=u'\u2212'
#mpl.rc('lines',linewidth=1.2)

sfctype_list=['180','181','182','183','187']
outputpath=rootpath+'/AlbanyWork/Prospectus/Experiments/HazyDA/Images/DiagFiles/conv/biasrms'
inputpath=rootarch

varlist=['sst'] #['ps','sst','gps','q','t','uv','tcp']
unitlist=['K'] #['mb','K','%','g/kg','K','m/s','mb']
bufrtype='all' # SST: 181-199
explist=np.array(['hazyda_ctrl','hazyda_aero'])
expnlist=['CTL','AER']
nexp=explist.size
exps_fname_str=''
for exp in expnlist:
    exps_fname_str += exp+'_'
exps_fname_str=exps_fname_str[:-1]

sensor = 'iasi_metop-a'
chkwvn = 962.5
sel_radqc = 13

sdate=2020061000
edate=2020071018
hint=6

loop='ges' #ges,anl
area='r2o10'
useqc=1

if (loop=='ges'):
   lpstr='OMF'
elif (loop=='anl'):
   lpstr='OMA'
if (useqc):
    qcflg='qc'
else:
    qcflg='noqc'
minlon, maxlon, minlat, maxlat, crosszero, cyclic=setarea.setarea(area)
print(area,minlat,maxlat,minlon,maxlon,crosszero)
cornerll = [minlat, maxlat, minlon, maxlon]

zpltlst=[0,1,2,3,4,5,6,7,8]

imgsavpath=outputpath+'/'+area
if ( not os.path.exists(imgsavpath) ):
   os.makedirs(imgsavpath)

date1 = pd.to_datetime(sdate,format='%Y%m%d%H')
date2 = pd.to_datetime(edate,format='%Y%m%d%H')
delta = timedelta(hours=6)
dates = pd.date_range(start=date1, end=date2, freq=delta)
tnum = dates.size

rule = rrulewrapper(DAILY, byhour=6, interval=5)
loc = RRuleLocator(rule)
formatter = DateFormatter('%Y %h %n %d %Hz')

print('Total cases number is %d' % tnum )

ptop=np.array((1000.,900.,800.,600.,400.,300.,250.,200.,150.,100.,50.,0.))
pbot=np.array((1200.,1000.,900.,800.,600.,400.,300.,250.,200.,150.,100.,50.))
pmid=0.5*(ptop+pbot)
znum=ptop.size
zlabels=[]
for z in np.arange(znum):
    zlabels.append('%.0i-%.0i'%(ptop[z],pbot[z]))

uidx=0
for var in varlist:
    unit=unitlist[uidx]
    sfcflag=(var=='ps' or var=='sst' or var=='tcp'
             or bufrtype in sfctype_list)
    if ( sfcflag ):
        icount=np.zeros((tnum,2))
    else:
        icount=np.zeros((tnum,znum,2))
    d=0
    for date in dates:
        datestr=date.strftime('%Y%m%d%H') 
        cnvdfile='diag_conv_'+var+'_'+loop+'.'+datestr+'.'+diagsuffix
        infile1=inputpath+'/'+explist[0]+'/'+datestr+'/'+cnvdfile
        infile2=inputpath+'/'+explist[1]+'/'+datestr+'/'+cnvdfile
        radfile=inputpath+'/'+explist[1]+'/'+datestr+'/diag_'+sensor+'_'+loop+'.'+str(date)+'.nc4'
        if (os.path.exists(infile1) and os.path.exists(infile2)):
            print('Processing Cnvfile: %s' %(cnvdfile))
            ds1=read_cnv_ncdiag(infile1,useqc=useqc,sel_bufr=bufrtype,is_sfc=sfcflag,area=area,cornerll=cornerll)
            ds2=read_cnv_ncdiag(infile2,useqc=useqc,sel_bufr=bufrtype,is_sfc=sfcflag,area=area,cornerll=cornerll)
            ds3=read_rad_ncdiag(radfile,chkwvn=chkwvn,sel_qc=sel_radqc,area=area,cornerll=cornerll)
            #print('Load Data Elapsed: %f [s]' %(end-start))
            try:
                omg_mean
            except NameError:
                rad_qccnt[:]=np.nan
                if (sfcflag):
                    omg_mean=np.zeros((tnum,2),dtype='float')
                    omg_rmsq=np.zeros_like(omg_mean)
                    omg_mean[:,:]=np.nan
                    omg_rmsq[:,:]=np.nan
                else:
                    omg_mean=np.zeros((tnum,znum,2),dtype='float')
                    omg_rmsq=np.zeros_like(omg_mean)
                    omg_mean[:,:,:]=np.nan
                    omg_rmsq[:,:,:]=np.nan
        else:
            print('%s is not existing'%(cnvdfile))
            d=d+1
            continue
        
        rad_qccnt[d]=ds3.omb.size
            
        if (sfcflag):
            dpar1=ds1.omb_nbc
            dpar2=ds2.omb_nbc
            
            icount[d,0]=dpar1.size
            icount[d,1]=dpar2.size
            
            omg_mean[d,0]=np.nanmean(dpar1)
            omg_mean[d,1]=np.nanmean(dpar2)
            
            omg_rmsq[d,0]=np.sqrt(np.nanmean(np.square(dpar1)))
            omg_rmsq[d,1]=np.sqrt(np.nanmean(np.square(dpar2)))
            
        else:
            dpar1=ds1.omb_nbc
            dpar2=ds2.omb_nbc
            
            pres1=ds1.pres
            pres2=ds2.pres
                
            for z in np.arange(znum):
                zmask1=((pres1<pbot[z])&(pres1>ptop[z]))
                zmask2=((pres2<pbot[z])&(pres2>ptop[z]))
                
                icount[d,z,0]=np.count_nonzero(zmask1)
                icount[d,z,1]=np.count_nonzero(zmask2)
                
                zdpar1=xa.where(zmask1,dpar1,np.nan)
                zdpar2=xa.where(zmask2,dpar2,np.nan)
                
                omg_mean[d,z,0]=np.nanmean(zdpar1)
                omg_mean[d,z,1]=np.nanmean(zdpar2)
            
                omg_rmsq[d,z,0]=np.sqrt(np.nanmean(np.square(zdpar1)))
                omg_rmsq[d,z,1]=np.sqrt(np.nanmean(np.square(zdpar2)))

        d=d+1

    aod_ts = opendap_m2_aod(sdate,edate,hint,area,cornerll=cornerll,varname='TOTEXTTAU',outfmt='ts')
    aod_ts_data = aod_ts.data
                
    if (sfcflag):
        fig,ax=plt.subplots(3,1,sharex=True,figsize=(9,5.7))
        fig.subplots_adjust(hspace=0.1)
        for a in np.arange(3):
            ax[a].set_prop_cycle(color=['blue','red','green'])
            ax[a].grid()
        
        ax[0].plot_date(dates,aod_ts_data,'k',lw=1)
        ax[0].set_ylabel('MERRA-2 AOD')
        ax[0].xaxis.set_major_locator(loc)
        ax[0].xaxis.set_major_formatter(formatter)
        ax[0].xaxis.set_tick_params(labelsize=10)
        ax[0].set_title('%s %s[%s]' %(area,var.upper(),unit),loc='left')
        ax[1].plot_date(dates,omg_rmsq,'--o',lw=1,ms=1.5)
        ax[1].set_ylabel('RMS %s [%s]' %(lpstr,unitlist[uidx]))
        ax[2].plot_date(dates,omg_mean,'-o',lw=1,ms=1.5)
        ax[2].set_ylabel('Mean %s [%s]'%(lpstr,unitlist[uidx]))
        lglist=np.zeros((2,2),dtype='<U30')
        for ex in np.arange(nexp):
            lglist[0,ex]=expnlist[ex]+'(%.2f)' %(np.nanmean(omg_rmsq[:,ex]))
            lglist[1,ex]=expnlist[ex]+'(%.2f)' %(np.nanmean(omg_mean[:,ex]))
        ax[1].legend(lglist[0,:])
        ax[2].legend(lglist[1,:])
        if (fsave):
            fname='%s/%s_%s_%s_%s_%s_bufr%s_BIASRMS_AOD.%s_%s.png'%(imgsavpath,area,loop,var,exps_fname_str,qcflg,bufrtype,sdate,edate)
            print(fname,flush=1)
            fig.savefig(fname, dpi=quality)
            plt.close()
        
    else:
        fig,ax=plt.subplots(1,3,sharey=True)
        fig.subplots_adjust(left=0.2,right=0.95,wspace=0.05)
        ax[0].set_prop_cycle(color=['blue','red'],linestyle=['-','-'],marker=['o','o'])
        ax[1].set_prop_cycle(color=['blue','red'],linestyle=['--','--'],marker=['o','o'])
        ax[2].set_prop_cycle(color=['k','k'],linestyle=['-','--'],marker=['o','o'])
        biasplot=np.nanmean(omg_mean,axis=0)
        rmsplot=np.nanmean(omg_rmsq,axis=0)
        diffplot=np.zeros_like(biasplot)
        diffplot[:,0]=np.diff(biasplot,axis=1)[:,0]
        diffplot[:,1]=np.diff(rmsplot,axis=1)[:,0]
        ax[0].plot(biasplot,pmid,lw=1,ms=1.5)
        ax[0].set_title('Mean %s [%s]'%(lpstr,unitlist[uidx]))
        ax[1].plot(rmsplot,pmid,lw=1,ms=1.5)
        ax[1].set_title('RMS %s [%s]'%(lpstr,unitlist[uidx]))
        ax[2].plot(diffplot,pmid,lw=1,ms=1.5)
        ax[2].set_title('%s%s%s [%s]'%(expnlist[1],minussign,expnlist[0],unit))
        ax[0].invert_yaxis()
        ax[0].set_yscale('log')
        ax[0].set_yticks(pmid[::2])
        ax[0].set_yticklabels(zlabels[::2])
        ax[0].set_ylabel('Pressure [hPa]')
        fig.suptitle('%s [%s]' %(var.upper(),unit))
        ax[0].legend(expnlist)
        ax[2].legend(['Mean','RMS'],loc=2)
        ax[0].grid()
        ax[1].grid()
        ax[2].grid()
        if (fsave):
            fname=('%s/%s_%s_%s_%s_%s_bufr%s_BIASRMS.%s_%s.png'
                   %(imgsavpath,area,loop,var,exps_fname_str,qcflg,bufrtype,sdate,edate))
            print(fname,flush=1)
            fig.savefig(fname, dpi=quality)
            plt.close()
        
        for z in zpltlst:
            fig,ax=plt.subplots(2,1,sharex=True,figsize=(9,3.8))
            fig.subplots_adjust(hspace=0.1)
            for a in np.arange(2):
               ax[a].set_prop_cycle(color=['blue','red'])
               ax[a].grid()
        
            ax[1].plot_date(dates,omg_mean[:,z,:],'-o',ms=1.5,lw=1)
            ax[0].xaxis.set_major_locator(loc)
            ax[0].xaxis.set_major_formatter(formatter)
            ax[0].xaxis.set_tick_params(labelsize=10)
            ax[0].set_title('%s %s [%s] %.1f %s %.1f [hPa]' %(area,var.upper(),unit,ptop[z],minussign,pbot[z]),loc='left')
            ax[0].plot_date(dates,omg_rmsq[:,z,:],'--o',ms=1.5,lw=1)
            ax[0].set_ylabel('RMS %s [%s]'%(lpstr,unitlist[uidx]))
            ax[1].set_ylabel('Mean %s [%s]'%(lpstr,unitlist[uidx]))
            lglist=np.zeros((2,2),dtype='<U30')
            for ex in np.arange(nexp):
               lglist[0,ex]=expnlist[ex]+'(%.2f)' %(np.nanmean(omg_rmsq[:,z,ex]))
               lglist[1,ex]=expnlist[ex]+'(%.2f)' %(np.nanmean(omg_mean[:,z,ex]))
            ax[0].legend(lglist[0,:])
            ax[1].legend(lglist[1,:])
            
            if (fsave):
               fname=('%s/%s_%s_%s_%s_%i_%s_bufr%s_BIASRMS_AOD.%s_%s.png' 
                      %(imgsavpath,area,loop,var,exps_fname_str,pbot[z],qcflg,bufrtype,sdate,edate))
               print(fname,flush=1)
               fig.savefig(fname, dpi=quality)
               plt.close()
                
    uidx=uidx+1
    

