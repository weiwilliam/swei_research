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

tlsize=12 ; lbsize=10
mpl.rc('axes', titlesize=tlsize,labelsize=lbsize)
mpl.rc('xtick',labelsize=lbsize)
mpl.rc('ytick',labelsize=lbsize)
mpl.rc('legend',fontsize='large')
fsave=1 ; ffmt='png' ; ptsize=4
axe_w=3 ; axe_h=3 ; quality=300

# Projection setting
proj=ccrs.PlateCarree(globe=None)

# Plotting setup
sdate=2020061000
edate=2020071018
hint=6
aertype='Dust'
exp='AerObserver'
sensor='iasi_metop-a'
spectral_range=slice(700,1300)
loop='ges' #ges,anl
#if loop=='anl':
#    tlstr='OMA'
#elif loop=='ges':
#    tlstr='OMF'
outtxt=0 # write dates output
plt2d=0  # plot single cycle 2d map
pltxy=0  # plot period statistics
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

# Data path setup
archpath=rootarch+'/Prospectus/AeroObsStats/nc_diag'
outpath=rootpath+'/AlbanyWork/Prospectus/Experiments/AeroObsStats/images/'
archdir=archpath+'/'+exp
savedir=outpath+'/'+exp+'/hist2d/'+aertype
if ( not os.path.exists(savedir) ):
    os.makedirs(savedir)

cat_lst=['aer_qc0','aer_qc7','aer_qcother']
clridx=[0,11,4,9]
cmap=setup_cmap('amwg',clridx)
aertypelvs=[0.5,1.5,2.5,3.5]
aer_norm = mpcrs.BoundaryNorm(aertypelvs,len(clridx)+1,extend='both')

chkwvn_list=[906.25]

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

# if (outtxt):
#     filename='%s/%s_%s.aer_qcflg_rate.%s.txt' %(savedir,str(sdate),str(edate),str(chkwvn))
#     datesfile=open(filename,'w')
#     datesfile.write('Date           QC=0          QC=7         Others \n')
    
dates_count=0
for date in dlist:
    raddfile='diag_'+sensor+'_'+loop+'.'+str(date)+'.nc4'
    infile1=archdir+'/'+str(date)+'/'+raddfile

    if (os.path.exists(infile1)):
        print('Processing Radfile: %s' %(raddfile))
        ds1=xa.open_dataset(infile1)
        npts=int(ds1.nobs.size/ds1.nchans.size)
        nchs=ds1.nchans.size
        ds1=ds1.assign_coords(nuchan=('wavenumber',ds1.wavenumber))
        ds1=ds1.swap_dims({"nchans":"wavenumber"}) #replace the dimension of channel by channel indices
        wavelength=1e+04/ds1.wavenumber
        chkwvn_list=ds1.wavenumber.sel(wavenumber=spectral_range)[ds1.use_flag.sel(wavenumber=spectral_range)==1]
        #usedchidx=np.where(ds1.iuse_rad==1)[0]
        #unusedchidx=np.where(ds1.iuse_rad==-1)[0]
        #wvldiff=abs(np.subtract(wavelength[usedchidx],chkwvl))
        #chkwvlidx=usedchidx[np.where(wvldiff==wvldiff.min())[0][0]]
        #print('Check wavelength: %.2f' %(wavelength[chkwvlidx]))
        dates_count+=1
    else:
        print('%s is not existing'%(raddfile))
        continue
    
    # Observation lat/lon
    rlat1=np.reshape(ds1.Latitude.values,(npts,nchs))
    rlon1=np.reshape(ds1.Longitude.values,(npts,nchs))
    qcflags=np.reshape(ds1.QC_Flag.values,(npts,nchs))
    obs1=np.reshape(ds1.Observation.values,(npts,nchs))
    sim1=np.reshape(ds1.Simulated_Tb.values,(npts,nchs))
    clr1=np.reshape(ds1.Clearsky_Tb.values,(npts,nchs))
    tmpds=xa.Dataset({'rlon1':(['obsloc'],rlon1[:,0]),
                      'rlat1':(['obsloc'],rlat1[:,0]),
                      'qcflag':(['obsloc','wavenumber'],qcflags),
                      'tb_obs':(['obsloc','wavenumber'],obs1),
                      'tb_sim':(['obsloc','wavenumber'],sim1),
                      'tb_clr':(['obsloc','wavenumber'],clr1)},
                     coords={'obsloc':np.arange(npts),
                             'wavenumber':ds1.wavenumber.values})
    
    tmpds=tmpds.sel(wavenumber=spectral_range)
    
    if (date==str(sdate)):
        ds_all=tmpds
    else:
        ds_all=xa.concat((ds_all,tmpds),dim='obsloc')


# if (outtxt):
#     datesfile.close()

# for chkwvn in [906.25]:
for chkwvn in chkwvn_list:
    ds_chk=ds_all.sel(wavenumber=chkwvn)
    tb_sim=ds_chk.tb_sim
    tb_clr=ds_chk.tb_clr
    tb_obs=ds_chk.tb_obs

    omb=tb_obs-tb_sim
    aereff_fg=tb_sim-tb_clr
    aereff_obs=tb_obs-tb_clr
    aereff=0.5*abs(aereff_fg)+0.5*abs(aereff_obs)
    
    qc0_msk=(ds_chk.qcflag==0.)
    qc7_msk=(ds_chk.qcflag==7.)
    qc13_msk=(ds_chk.qcflag==13.)

    pltmask=(qc0_msk)|(qc13_msk)
    if (plthist2d):
        for pltset in [1,2]:
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
                y_label='O-B [K]'
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
                
            if (pltset==2):
                bndry_line_ny=np.arange(-30,-2)
                bndry_line_nx=bndry_line_ny/1.8
                bndry_line_py=np.arange(3,31)
                bndry_line_px=bndry_line_py/1.8
                plt.plot(abs(bndry_line_nx),bndry_line_ny,'k',linestyle='dashed')
                plt.plot(bndry_line_px,bndry_line_py,'k',linestyle='dashed')
                
            ax.set_xlabel(x_label)
            ax.set_ylabel(y_label)
            
            if (fsave):
                fname=('2DPDF_%s_%s_%.2f.set%s.%s'
                        %(area,sensor,chkwvn,pltset,ffmt))
                fig.savefig(savedir+'/'+fname,dpi=quality)
                plt.close()
