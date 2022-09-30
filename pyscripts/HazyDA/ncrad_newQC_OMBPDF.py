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
sys.path.append(rootgit+'/pyscripts/functions')
import numpy as np
import xarray as xa
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mpcrs
import cartopy.crs as ccrs
import setuparea as setarea
from plot_utils import setupax_2dmap, plt_x2y, set_size
from utils import ndate,setup_cmap
from datetime import datetime, timedelta
import scipy.stats
from gsi_ncdiag import read_rad_ncdiag,read_rad_ncdiag0

tlsize=12 ; lbsize=10
mpl.rc('axes', titlesize=tlsize,labelsize=lbsize)
mpl.rc('xtick',labelsize=lbsize)
mpl.rc('ytick',labelsize=lbsize)
mpl.rc('legend',fontsize='small')
fsave=1 ; ffmt='png' ; ptsize=4
axe_w=3 ; axe_h=3 ; quality=300

# Plotting setup
sdate=2020061000
edate=2020092118
aertype='All'
hint=6
explist=['ncepgdas','observer_v2qc']
expnlist=['GDAS','w/o CldQC']
saved_suffix='v1qc'
sensor='iasi_metop-a'
chkwvn_list=[906.25,943.25,962.5,1096.0]
nchs=616
lutver='v5'
loop='ges' #ges,anl
usebc=0
pltpdf=1 # plot pdf

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
lutpath=rootpath+'/AlbanyWork/Prospectus/Experiments/AeroObsStats/SD_LUT'
outpath=rootpath+'/AlbanyWork/Prospectus/Experiments/HazyDA/Images/DiagFiles/rad/observer'
archdir0=rootarch+'/'+explist[0]
archdir1=rootarch+'/'+explist[1]

savedir=outpath+'/pdf_'+saved_suffix
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
    raddfile0='diag_'+sensor+'_'+loop+'.'+str(date)+'.nc'
    raddfile1='diag_'+sensor+'_'+loop+'.'+str(date)+'.nc4'
    infile0=archdir0+'/'+str(date)+'/'+raddfile0
    infile1=archdir1+'/'+str(date)+'/'+raddfile1

    if (os.path.exists(infile0) and 
        os.path.exists(infile1)):
        print('Processing Radfile: %s' %(raddfile1),flush=1)
        dates_count+=1
    else:
        print('%s is not existing'%(raddfile1),flush=1)
        continue

    tmpds0=read_rad_ncdiag0(infile0,chkwvn=chkwvn_list,cal_ae=1,get_water_frac=1)
    tmpds1=read_rad_ncdiag(infile1, chkwvn=chkwvn_list,cal_ae=1,get_water_frac=1)
    
    if (date==str(sdate)):
        ds_all0=tmpds0
        ds_all1=tmpds1
    else:
        ds_all0=xa.concat((ds_all0,tmpds0),dim='obsloc')
        ds_all1=xa.concat((ds_all1,tmpds1),dim='obsloc')

total_obscounts=ds_all1.obsloc.size
ds_all1=ds_all1.assign_coords(obsloc=np.arange(total_obscounts))

#satinfo_csv=lutpath+'/'+sensor+'_'+str(nchs)+'_stats.'+lutver+'.csv'
#lutdf=pd.read_csv(satinfo_csv)
#filter = ((lutdf.Aer_sen==1.)&(lutdf.iuse==1.)&
#          (lutdf.wavenumber>=spectral_range.start)&
#          (lutdf.wavenumber<=spectral_range.stop))
#tmpdf=lutdf.loc[filter,:]
#chkwvn_list=tmpdf.wavenumber.values

binsize=0.1
halfbin=0.5*binsize
hist_x_edge=np.arange(-6.-halfbin,6.+binsize,binsize)
bin_center=(hist_x_edge+halfbin)[:-1]

for chkwvn in chkwvn_list:
    ds_chk0=ds_all0.sel(wavenumber=chkwvn)
    tb_obs0=ds_chk0.tb_obs
    # varinv1=ds_chk1.varinv
    
    ds_chk1=ds_all1.sel(wavenumber=chkwvn)
    tb_obs1=ds_chk1.tb_obs
    tb_sim1=ds_chk1.tb_sim
    tb_clr1=ds_chk1.tb_clr
    #varinv1=ds_chk1.varinv
    
    if (usebc):
        omb0=ds_chk0.omb_bc
        omb1=ds_chk1.omb_bc
        omb_clr1=(tb_obs1-tb_clr1)-(omb1-tb_obs1+tb_sim1)
    else:
        omb0=ds_chk0.omb_nbc
        omb1=ds_chk1.omb_nbc
        omb_clr1=tb_obs1-tb_clr1
    
    water_msk0=(ds_chk0.water_frac>0.999)
    good_msk0=(ds_chk0.qcflag==0.)
    gross_msk0=(ds_chk0.qcflag==3.)
    cld_msk0=(ds_chk0.qcflag==7.)
    tzr_msk0=(ds_chk0.qcflag==10.)
    sfcir_msk0=(ds_chk0.qcflag==53.)
    passed_msk0=(good_msk0)
    # ori_msk=(good_msk)|(aer_msk)|(bust_msk)|(tzr_msk)

    water_msk1=(ds_chk1.water_frac>0.999)
    good_msk1=(ds_chk1.qcflag==0.)
    gross_msk1=(ds_chk1.qcflag==3.)
    cld_msk1=(ds_chk1.qcflag==7.)
    tzr_msk1=(ds_chk1.qcflag==10.)
    sfcir_msk1=(ds_chk1.qcflag==53.)
    bust_msk1=(ds_chk1.qcflag==55.)
    aercld_msk1=(ds_chk1.qcflag==57.)
    if (saved_suffix=='v1qc'):
       bustqc=(abs(omb1)>3.)&(abs(omb1)>1.8*ds_chk1.Ae)
       tmpmsk1=(aercld_msk1)&(~bustqc)
       omb1=xa.where(tmpmsk1,omb_clr1,omb1)
       aer_msk1=(ds_chk1.qcflag==13.)|(tmpmsk1)
    else:
       aer_msk1=(ds_chk1.qcflag==13.)
    passed_msk1=(good_msk1)|(aer_msk1)
    # ori_msk=(good_msk)|(aer_msk)|(bust_msk)|(tzr_msk)
    
    pltmsk0=passed_msk0 ; cnts0=np.count_nonzero(pltmsk0)
    pltmsk1=passed_msk1 ; cnts1=np.count_nonzero(pltmsk1)

    cntslist=[cnts0,cnts1]
    leglist=[]
    for lg in zip(expnlist,cntslist):
        cntstr="{:,}".format(lg[1])
        lgstr='%s (%s)' %(lg[0],cntstr)
        leglist.append(lgstr)

    if (pltpdf):            
        pltda_x0=omb0[pltmsk0==1]
        pltda_x1=omb1[pltmsk1==1]
        
        x_label='OMB [K]'
        hdata0, tmpbins=np.histogram(pltda_x0, bins=hist_x_edge, density=0)
        hdata1, tmpbins=np.histogram(pltda_x1, bins=hist_x_edge, density=0)
    
        omb_mean=np.zeros_like(bin_center,dtype='float')
        omb_sd=np.zeros_like(bin_center,dtype='float')
        # counts=np.zeros_like(bin_center,dtype='int')
        
        tistr=('%s (%.2f $cm^{-1}$)' %(sensor,chkwvn))
        
        fig=plt.figure()
        ax=plt.subplot()
        set_size(axe_w,axe_h,l=0.18,r=0.88)
        ax.set_title(tistr,loc='left')
        ax.set_xlabel(x_label)
        ax.plot(bin_center,hdata0,'tab:blue')
        ax.plot(bin_center,hdata1,'tab:red')
        ax.set_ylabel("Counts")
        ax.vlines(0.,0,1,transform=ax.get_xaxis_transform(),colors='grey',linestyle='dashed',linewidth=0.7)
        ax.legend(leglist,loc=2)
        
        if (fsave):
            fname=('%s/PDF_OMB_%s_%s_%s.%.2f.%s'
                    %(savedir,area,sensor,bcflg,chkwvn,ffmt))
            print(fname,flush=1)
            fig.savefig(fname,dpi=quality)
            plt.close()
