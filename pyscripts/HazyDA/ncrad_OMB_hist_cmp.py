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
    rootpath='/data/users/swei/'
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
import scipy.stats
from gsi_ncdiag import read_rad_ncdiag

tlsize=12 ; lbsize=10
mpl.rc('axes', titlesize=tlsize,labelsize=lbsize)
mpl.rc('xtick',labelsize=lbsize)
mpl.rc('ytick',labelsize=lbsize)
mpl.rc('legend',fontsize='small')
fsave=1 ; ffmt='png' ; ptsize=4
axe_w=3 ; axe_h=3 ; quality=300

# Plotting setup
sdate=2020061000
edate=2020071018
hint=6
exp='hazyda_aero'
expname='AER'
sensor='iasi_metop-a'
spectral_range=slice(600,1300)
chkwvn_list=[962.5]
loop='ges' #ges,anl
if loop=='anl':
    tlstr='OMA'
elif loop=='ges':
    tlstr='OMF'
plthist=1 # plot 2d histogram

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


# Data path setup
outpath=rootpath+'/AlbanyWork/Prospectus/Experiments/HazyDA/Images'
archdir0=rootarch+'/'+exp
savedir=outpath+'/DiagFiles/rad/pdf/'+area
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
    raddfile0='diag_'+sensor+'_'+loop+'.'+str(date)+'.nc4'
    infile0=archdir0+'/'+str(date)+'/'+raddfile0
    infile1=archdir1+'/'+str(date)+'/'+raddfile1

    if ( os.path.exists(infile0) ): 
        print('Processing Radfile: %s' %(raddfile0),flush=1)
        dates_count+=1
    else:
        print('%s is not existing'%(raddfile0),flush=1)
        continue
    
    # Observation lat/lon from exp 0 (baseline)
    tmpds0=read_rad_ncdiag(chkwvn=chkwvn_list)
    tmpds1=read_rad_ncdiag(chkwvn=chkwvn_list)
    
    if (date==str(sdate)):
        ds_all0=tmpds0
        ds_all1=tmpds0
    else:
        ds_all0=xa.concat((ds_all0,tmpds0),dim='obsloc')
        ds_all1=xa.concat((ds_all1,tmpds0),dim='obsloc')

total_obscounts=ds_all0.obsloc.size
ds_all0=ds_all0.assign_coords(obsloc=np.arange(total_obscounts))
ds_all1=ds_all1.assign_coords(obsloc=np.arange(total_obscounts))

binsize=0.1
halfbin=0.5*binsize
hist_x_edge=np.arange(-10,10.+binsize,binsize)
bin_center=(hist_x_edge+halfbin)[:-1]

if type(chkwvn_list)==slice:
   chkwvn_list=ds_all0.wavenumber

for chkwvn in chkwvn_list:
    ds_chk0=ds_all0.sel(wavenumber=chkwvn)
    ds_chk1=ds_all1.sel(wavenumber=chkwvn)

    omb_bc0=ds_chk0.omb_bc
    omb_nbc0=ds_chk0.omb_nbc
    omb_bc1=ds_chk1.omb_bc
    omb_nbc1=ds_chk1.omb_nbc
    
    good_msk0=(ds_chk0.qcflag==0.)
    gross_msk0=(ds_chk0.qcflag==3.)
    cld_msk0=(ds_chk0.qcflag==7.)
    tzr_msk0=(ds_chk0.qcflag==10.)
    aer_msk0=(ds_chk0.qcflag==13.)
    sfcir_msk0=(ds_chk0.qcflag==53.)
    bust_msk0=(ds_chk0.qcflag==55.)
    aercld_msk0=(ds_chk0.qcflag==57.)

    good_msk1=(ds_chk1.qcflag==0.)
    gross_msk1=(ds_chk1.qcflag==3.)
    cld_msk1=(ds_chk1.qcflag==7.)
    tzr_msk1=(ds_chk1.qcflag==10.)
    aer_msk1=(ds_chk1.qcflag==13.)
    sfcir_msk1=(ds_chk1.qcflag==53.)
    bust_msk1=(ds_chk1.qcflag==55.)
    aercld_msk1=(ds_chk1.qcflag==57.)

    passed_msk0=(good_msk0)
    passed_msk1=(good_msk1)|(aer_msk1)
    
    if (plthist):            
        bcflg1='bc'
        bcflg2='nobc'
        x_label='OMB'
        y_label='Data Count'
        pltda_x1=omb_bc0[good_msk0]
        pltda_x2=omb_bc0[aer_msk0]
        hdata1, tmpbins=np.histogram(pltda_x1, bins=hist_x_edge, density=0)
        hdata2, tmpbins=np.histogram(pltda_x2, bins=hist_x_edge, density=0)
        pltda_x3=omb_nbc0[good_msk0]
        pltda_x4=omb_nbc0[aer_msk0]
        hdata3, tmpbins=np.histogram(pltda_x3, bins=hist_x_edge, density=0)
        hdata4, tmpbins=np.histogram(pltda_x4, bins=hist_x_edge, density=0)
        tmpymax=np.nanmax((hdata1,hdata2,hdata3,hdata4)) 
        ymax=np.ceil(tmpymax/np.power(10,int(np.log10(tmpymax))))*np.power(10,int(np.log10(tmpymax)))
        #print('ymax at %.2f = %s, %s' %(chkwvn,tmpymax,ymax) )

        #omb_mean=np.zeros_like(bin_center,dtype='float')
        #omb_sd=np.zeros_like(bin_center,dtype='float')
        # counts=np.zeros_like(bin_center,dtype='int')
        
        tistr=('%s (%.2f $\mathrm{cm^{-1}}$)' %(sensor,chkwvn))
        
        fig=plt.figure()
        ax=plt.subplot()
        set_size(axe_w,axe_h,l=0.2,r=0.9)
        ax.set_title(tistr,loc='left')
        ax.set_xlabel(x_label)
        ax.plot(bin_center,hdata1,'bo',fillstyle='none',ms=2.5)
        ax.plot(bin_center,hdata2,'ro',fillstyle='none',ms=2.5)
        #ax.set_yscale("log")
        #ax.set_ylabel("PDF")
        ax.set_ylim(-int(ymax*0.05),ymax)
        ax.set_ylabel(y_label)
        ax.vlines(0.,0,1,transform=ax.get_xaxis_transform(),colors='grey',linestyle='dashed',linewidth=0.7)
        ax.legend(['Clear','Hazy'],loc=2)
        
        if (fsave):
            fname=('%s/PDF_OMB_%s_%s_%s_%.2f_%s.%s_%s.%s'
                    %(savedir,area,expname,sensor,chkwvn,bcflg1,sdate,edate,ffmt))
            print(fname,flush=1)
            fig.savefig(fname,dpi=quality)
            plt.close()

        fig=plt.figure()
        ax=plt.subplot()
        set_size(axe_w,axe_h,l=0.2,r=0.9)
        ax.set_title(tistr,loc='left')
        ax.set_xlabel(x_label)
        ax.plot(bin_center,hdata3,'bo',fillstyle='none',ms=2.5)
        ax.plot(bin_center,hdata4,'ro',fillstyle='none',ms=2.5)
        #ax.set_yscale("log")
        #ax.set_ylabel("PDF")
        ax.set_ylim(-int(ymax*0.05),ymax)
        ax.set_ylabel(y_label)
        ax.vlines(0.,0,1,transform=ax.get_xaxis_transform(),colors='grey',linestyle='dashed',linewidth=0.7)
        ax.legend(['Clear','Hazy'],loc=2)
        
        if (fsave):
            fname=('%s/PDF_OMB_%s_%s_%s_%.2f_%s.%s_%s.%s'
                    %(savedir,area,expname,sensor,chkwvn,bcflg2,sdate,edate,ffmt))
            print(fname,flush=1)
            fig.savefig(fname,dpi=quality)
            plt.close()
