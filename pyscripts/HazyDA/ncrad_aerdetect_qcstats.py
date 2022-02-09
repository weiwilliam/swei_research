# -*- coding: utf-8 -*-
"""
Created on Tue Jun 29 16:06:23 2021

@author: ck102

Aerosol detection based on CADS 3.1 from NWP SAF

"""
import sys, os, platform
import numpy as np
import xarray as xa
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mpcrs
import cartopy.crs as ccrs
os_name=platform.system()
if (os_name=='Darwin'):
    rootpath='/Users/weiwilliam'
    rootarch='/Volumes/WD2TB/ResearchData'
elif (os_name=='Windows'):
    rootpath='F:\GoogleDrive_NCU\Albany'
    rootarch='F:\ResearchData'
sys.path.append(rootpath+'/AlbanyWork/Utility/Python3/functions')
import setuparea as setarea
from plot_utils import setupax_2dmap, plt_x2y, set_size
from utils import ndate,setup_cmap
from datetime import datetime, timedelta

tlsize=12 ; lbsize=12
mpl.rc('axes', titlesize=tlsize,labelsize=lbsize)
mpl.rc('xtick',labelsize=lbsize)
mpl.rc('ytick',labelsize=lbsize)
mpl.rc('legend',fontsize='xx-large')
fsave=1 ; ffmt='png' ; ptsize=4
axe_w=6 ; axe_h=3 ; quality=300

# Projection setting
proj=ccrs.PlateCarree(globe=None)
"""
&Aerosol_Detect_Coeffs
N__Num_Aerosol_Tests = 3
N__Num_Aerosol_Chans(1:3) = 4, 4, 4
N__Aerosol_Chans(1,1:4) = 1340, 2348, 1782, 2356
N__Aerosol_Chans(2,1:4) = 2093, 2348, 1782, 2356
N__Aerosol_Chans(3,1:4) = 1782, 2356, 1782,  753
N__Mean_Aerosol_Chans = 11
R__Aerosol_TBD(1,1:2) = 0.2, -1.55
R__Aerosol_TBD(2,1:2) = -1.5, 999.9
R__Aerosol_TBD(3,1:2) = -1.8, -2.0
R__Rank_Thres_Coeff(1:3) = -0.01, 2.1, -3.9
R__Unclassified_Thres = 0.4
R__Land_Fraction_Thres = 0.5
"""
# Plotting setup
sdate=2020082200
edate=2020092118
hint=6
exp='GDAS'
sensor='iasi_metop-a'
chkwvn=906.25
loop='ges' #ges,anl
#if loop=='anl':
#    tlstr='OMA'
#elif loop=='ges':
#    tlstr='OMF'
outtxt=1 # write dates output
plt2d=1  # plot single cycle 2d map
pltxy=0  # plot period statistics

spec_int=0.25

area='Glb'
minlon, maxlon, minlat, maxlat, crosszero, cyclic=setarea.setarea(area)
print(minlat,maxlat,minlon,maxlon,crosszero,cyclic)
if (area=='Glb'):
   minlon=-180. ; maxlon=180.
cornll=[minlat,maxlat,minlon,maxlon]

cbori='horizontal' #vertical, horizontal
if (cbori=='vertical'):
   cb_frac=0.025
   cb_pad=0.06
elif (cbori=='horizontal'):
   cb_frac=0.04
   cb_pad=0.1

# Data path setup
archpath=rootarch+'/Prospectus/AeroObsStats/nc_diag'
outpath=rootpath+'/AlbanyWork/Prospectus/Experiments/AeroObsStats/images/'
archdir=archpath+'/'+exp
savedir=outpath+'/2dmap/aer_qcstat/'+exp
if ( not os.path.exists(savedir) ):
    os.makedirs(savedir)

cat_lst=['aer_qc0','aer_qc7','aer_qcother']
clridx=[0,11,4,9]
cmap=setup_cmap('amwg',clridx)
aertypelvs=[0.5,1.5,2.5,3.5]
aer_norm = mpcrs.BoundaryNorm(aertypelvs,len(clridx)+1,extend='both')

# Key channels
# 980 cm-1, 1232 cm-1, 1090.5 cm-1, 1234 cm-1, 1168 cm-1, and 833 cm-1,
# Setup for sensor
num_use_chans=4
ntest_wvn=np.zeros((num_use_chans),dtype='int')
thres_btd=np.zeros((2),dtype='float')
if ('iasi' in sensor):
    ntest_wvn[:4]=[980, 1232, 1090.5, 1234]
    nmean_chans=11
    thres_btd[:2]=[0.2, -1.55]

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

if (outtxt):
    filename='%s/%s_%s.aer_qcflg_rate.%s.txt' %(savedir,str(sdate),str(edate),str(chkwvn))
    datesfile=open(filename,'w')
    datesfile.write('Date           QC=0          QC=7         Others \n')
    
dates_count=0
for date in dlist:
    raddfile='diag_'+sensor+'_'+loop+'.'+str(date)+'.nc'
    infile1=archdir+'/'+raddfile

    if (os.path.exists(infile1)):
        print('Processing Radfile: %s' %(raddfile))
        ds1=xa.open_dataset(infile1)
        ds1=ds1.assign_coords(nuchan=('wavenumber',ds1.wavenumber))
        ds1=ds1.swap_dims({"channel":"wavenumber"}) #replace the dimension of channel by channel indices
        npts=ds1.obsloc.size  #int(ds1.nobs.size/ds1.nchans.size)
        nchs=ds1.channel.size #nchans.size
        wavelength=1e+04/ds1.wavenumber
        usedchidx=np.where(ds1.iuse_rad==1)[0]
        unusedchidx=np.where(ds1.iuse_rad==-1)[0]
        #wvldiff=abs(np.subtract(wavelength[usedchidx],chkwvl))
        #chkwvlidx=usedchidx[np.where(wvldiff==wvldiff.min())[0][0]]
        #print('Check wavelength: %.2f' %(wavelength[chkwvlidx]))
        dates_count+=1
    else:
        print('%s is not existing'%(raddfile))
        continue
    
    # Observation lat/lon
    rlat1=ds1.locinfo[0,:]
    rlon1=ds1.locinfo[1,:]
    
    # QC flag in IR window region
    qc0_tmp=(ds1.qcflag.sel(wavenumber=chkwvn)==0.)
    qc7_tmp=(ds1.qcflag.sel(wavenumber=chkwvn)==7.)
    
    # aerosol detection
    for i in np.arange(num_use_chans):
        keywvn=ntest_wvn[i]
        lb_wvn=keywvn-int(nmean_chans/2)*spec_int
        ub_wvn=keywvn+int(nmean_chans/2)*spec_int
        tmp_tb=ds1.tb_obs.sel(wavenumber=slice(lb_wvn,ub_wvn))
        if (date==str(sdate)):
            print('key wavenumber: %.2f, range: %.2f-%.2f, avail ch.: %i' 
                  %(keywvn,lb_wvn,ub_wvn,tmp_tb.wavenumber.size))
            print(tmp_tb.wavenumber)
        if (i==0):
            mean_bt=tmp_tb.mean(dim='wavenumber')
        else:
            mean_bt=xa.concat((mean_bt,tmp_tb.mean(dim='wavenumber')),dim='key_chan')

    btd1_tmp=mean_bt[0,:]-mean_bt[1,:]
    btd2_tmp=mean_bt[2,:]-mean_bt[3,:]
    aer_tmp=(btd1_tmp<=thres_btd[0])&(btd2_tmp<=thres_btd[1])
    
    aer_qc0=(aer_tmp)&(qc0_tmp)
    aer_qc7=(aer_tmp)&(qc7_tmp)
    aer_qcelse=(aer_tmp)&(~(qc0_tmp|qc7_tmp))

    if (plt2d):
        categories=xa.zeros_like(rlat1,dtype='float')
        categories[aer_qc0]=1.
        categories[aer_qc7]=2.
        categories[aer_qcelse]=3.
            
        fig,ax=setupax_2dmap(cornll,area,proj,lbsize=lbsize)
        set_size(axe_w,axe_h,l=0.1,r=0.85)
        pltmsk=(aer_tmp)
        #plttotal=np.count_nonzero(pltmsk)
        sc=ax.scatter(rlon1[pltmsk==1],rlat1[pltmsk==1],c=categories[pltmsk==1],
                      s=ptsize,cmap=cmap,norm=aer_norm,alpha=0.8,edgecolors='None')
        cidx=1
        cblblst=[]
        for cat in cat_lst:
            catcnts=np.count_nonzero(categories==cidx)
            cblblst.append('%s[%i]'%(cat,catcnts))
            cidx+=1
            
        cb=plt.colorbar(sc,orientation=cbori,fraction=cb_frac,ticks=[1.,2.,3.],pad=cb_pad,aspect=40)
        cb.ax.set_xticklabels(cblblst)
        #ax.set_title('%s (%i)' %(qcstr,np.count_nonzero(aer_tmp)),loc='left')
        if (fsave):
            fname='%s_%s_%s_%.2f_aer_qcflg_%s.%s'%(sensor,area,exp,chkwvn,date,ffmt)
            fig.savefig(savedir+'/'+fname,dpi=quality)
            plt.close()
    
    if (outtxt):
        aer_in_qc0=np.count_nonzero(aer_qc0)/np.count_nonzero(qc0_tmp)*100.
        aer_in_qc7=np.count_nonzero(aer_qc7)/np.count_nonzero(qc7_tmp)*100.
        aer_in_qcelse=np.count_nonzero(aer_qcelse)/np.count_nonzero(~(qc0_tmp|qc7_tmp))*100.
        outstr='%s %i %i %.2f%% %i %i %.2f%% %i %i %.2f%%\n' %(date,np.count_nonzero(aer_qc0),np.count_nonzero(qc0_tmp),aer_in_qc0,
                                                  np.count_nonzero(aer_qc7),np.count_nonzero(qc7_tmp),aer_in_qc7,
                                                  np.count_nonzero(aer_qcelse),np.count_nonzero(~(qc0_tmp|qc7_tmp)),
                                                  aer_in_qcelse)
        datesfile.write(outstr)
        
    if (date==str(sdate)):
        btd1_period=btd1_tmp
        btd2_period=btd2_tmp
        aer_msk=aer_tmp
        qc0_msk=qc0_tmp
        qc7_msk=qc7_tmp
    else:
        btd1_period=xa.concat((btd1_period,btd1_tmp),dim='obsloc')
        btd2_period=xa.concat((btd2_period,btd2_tmp),dim='obsloc')
        aer_msk=xa.concat((aer_msk,aer_tmp),dim='obsloc')
        qc0_msk=xa.concat((qc0_msk,qc0_tmp),dim='obsloc')
        qc7_msk=xa.concat((qc7_msk,qc7_tmp),dim='obsloc')
    
    del(mean_bt)

datesfile.close()

#aer_msk=(btd1_period<=thres_btd[0])&(btd2_period<=thres_btd[1])
noaer_msk=~aer_msk
print('aerosol counts:%i' %(np.count_nonzero(aer_msk)))
print('no aerosol counts:%i' %(np.count_nonzero(noaer_msk)))
print('Total counts: %i' %(aer_msk.size))

#qc_aer_msk=(qc_msk)&(aer_msk)
#print('aerosol-affected %s at %.2f(um): %i' 
#      %(qcstr,chkwvn,np.count_nonzero(qc_aer_msk.sel(wavenumber=chkwvn))))

#qc_counts_ave=np.zeros(qc_msk.wavenumber.size,dtype='float')
#for i in np.arange(qc_msk.wavenumber.size):
#    qc_counts_ave[i]=np.count_nonzero(qc_msk[i,:])/dates_count

