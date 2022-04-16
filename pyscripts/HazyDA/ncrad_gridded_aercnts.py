# -*- coding: utf-8 -*-
import sys, os, platform
machine='S4'
if (machine=='MBP'):
    rootpath='/Users/weiwilliam'
    rootarch='/Volumes/WD2TB/ResearchData'
elif (machine=='Desktop'):
    rootpath='F:\GoogleDrive_NCU\Albany'
    rootarch='F:\ResearchData'
    rootgit='F:\GitHub\swei_research'
elif (machine=='Hera'):
    rootpath='/scratch2/BMC/gsd-fv3-dev/Shih-wei.Wei'
    rootarch='/scratch2/BMC/gsd-fv3-dev/Shih-wei.Wei/ResearchData'
    rootgit='/home/Shih-wei.Wei/research'
elif (machine=='Cheyenne'):
    rootpath='/glade/work/swei/output/images'
    rootarch='/scratch2/BMC/gsd-fv3-dev/Shih-wei.Wei/ResearchData'
    rootgit='/glade/u/home/swei/research'
elif (machine=='S4'):
    rootpath='/data/users/swei'
    rootarch='/scratch/users/swei/ncdiag'
    rootgit='/home/swei/research'
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
from matplotlib.dates import (DAILY, DateFormatter,
                              rrulewrapper, RRuleLocator)

tlsize=12 ; lbsize=12
mpl.rc('axes', titlesize=tlsize,labelsize=lbsize)
mpl.rc('xtick',labelsize=lbsize)
mpl.rc('ytick',labelsize=lbsize)
mpl.rc('legend',fontsize='large')
fsave=1 ; ffmt='png' ; ptsize=4
axe_w=6 ; axe_h=3 ; quality=300

# Projection setting
proj=ccrs.PlateCarree(globe=None)

outputpath=rootpath+'/archive/HazyDA/gridded_diag'
inputpath=rootarch

# Plotting setup
expname='hazyda_ctrl'
sensor='iasi_metop-a'
sdate=2020082200
edate=2020092118
hint=6
sensor='iasi_metop-a'
selwvn=962.5
loop='ges' #ges,anl
degres=2.5
tkfreq=1
gen_data=0
gen_file=0
gen_plot=1

area='Glb'
minlon, maxlon, minlat, maxlat, crosszero, cyclic=setarea.setarea(area)
print(minlat,maxlat,minlon,maxlon,crosszero,cyclic)
if (area=='Glb'):
   minlon=-180. ; maxlon=180.
else:
   minlon=(minlon+180)%360-180
   maxlon=(maxlon+180)%360-180
cornll=[minlat,maxlat,minlon,maxlon]

cbori='vertical' #vertical, horizontal
if (cbori=='vertical'):
   cb_frac=0.025
   cb_pad=0.06
elif (cbori=='horizontal'):
   cb_frac=0.04
   cb_pad=0.1

syy=int(str(sdate)[:4]); smm=int(str(sdate)[4:6])
sdd=int(str(sdate)[6:8]); shh=int(str(sdate)[8:10])
eyy=int(str(edate)[:4]); emm=int(str(edate)[4:6])
edd=int(str(edate)[6:8]); ehh=int(str(edate)[8:10])

date1 = datetime(syy,smm,sdd,shh)
date2 = datetime(eyy,emm,edd,ehh)
delta = timedelta(hours=hint)
dates = pd.date_range(start=date1, end=date2, freq=delta)

rule = rrulewrapper(DAILY, byhour=6, interval=4)
loc = RRuleLocator(rule)
formatter = DateFormatter('%Y %h %n %d %Hz')

if (gen_data):
# Key channels
# 980 cm-1, 1232 cm-1, 1090.5 cm-1, 1234 cm-1, 1168 cm-1, and 833 cm-1,
# Setup for sensor
    num_use_chans=4
    ntest_wvn=np.zeros((num_use_chans),dtype='int')
    thres_btd=np.zeros((2),dtype='float')
    if ('iasi' in sensor):
        spec_int=0.25
        ntest_wvn[:4]=[980, 1232, 1090.5, 1234]
        nmean_chans=11
        thres_btd[:2]=[0.2, -1.55]
    
    tnum=0
    dlist=[]
    cdate=sdate
    while (cdate<=edate):
        dlist.append(str(cdate))
        tnum=tnum+1
        cdate=ndate(hint,cdate)
    
    didx=0
    for date in dlist:
        didx+=1
        cur_date=dates[didx-1]
        raddfile='diag_'+sensor+'_'+loop+'.'+str(date)+'.nc4'
        infile1=inputpath+'/'+expname+'/'+str(date)+'/'+raddfile
    
        if (os.path.exists(infile1)):
            print('Processing Radfile: %s' %(raddfile))
            ds1=xa.open_dataset(infile1)
            npts=int(ds1.nobs.size/ds1.nchans.size)
            nchs=ds1.nchans.size
            ds1=ds1.swap_dims({"nchans":"wavenumber"}) #replace the dimension of channel by channel indices
            
            rlat1=np.reshape(ds1.Latitude.values,(npts,nchs))
            rlon1=np.reshape(ds1.Longitude.values,(npts,nchs))
            rlon1=(rlon1+180)%360-180
            qcflags=np.reshape(ds1.QC_Flag.values,(npts,nchs))
            sim1=np.reshape(ds1.Simulated_Tb.values,(npts,nchs))
            sim_nbc1=np.reshape(ds1.Obs_Minus_Forecast_unadjusted.values,(npts,nchs))
            obs1=sim_nbc1+sim1
            tmpds=xa.Dataset({'rlon':(['obsloc'],rlon1[:,0]),
                              'rlat':(['obsloc'],rlat1[:,0]),
                              'qcflag':(['obsloc','wavenumber'],qcflags),
                              'tb_obs':(['obsloc','wavenumber'],obs1)},
                             coords={'obsloc':np.arange(npts),
                                     'wavenumber':ds1.wavenumber.values})
            # aerosol detection
            for i in np.arange(num_use_chans):
                keywvn=ntest_wvn[i]
                lb_wvn=keywvn-int(nmean_chans/2)*spec_int
                ub_wvn=keywvn+int(nmean_chans/2)*spec_int
                tmp_tb=tmpds.tb_obs.sel(wavenumber=slice(lb_wvn,ub_wvn))
                #if (date==str(sdate)):
                #    print('key wavenumber: %.2f, range: %.2f-%.2f, avail ch.: %i'
                #          %(keywvn,lb_wvn,ub_wvn,tmp_tb.wavenumber.size))
                #    print(tmp_tb.wavenumber)
                if (i==0):
                    mean_bt=tmp_tb.mean(dim='wavenumber')
                else:
                    mean_bt=xa.concat((mean_bt,tmp_tb.mean(dim='wavenumber')),dim='key_chan')
        
            btd1_tmp=mean_bt[0,:]-mean_bt[1,:]
            btd2_tmp=mean_bt[2,:]-mean_bt[3,:]
            aer_tmp=(btd1_tmp<=thres_btd[0])&(btd2_tmp<=thres_btd[1])
        
            tmpds=tmpds.sel(wavenumber=selwvn)
            aermsk_tmpds=xa.Dataset({'rlon':(['obsloc'],tmpds.rlon.data),
                                     'rlat':(['obsloc'],tmpds.rlat.data),
                                     'aer_msk':(['obsloc'],aer_tmp.data)},
                                    coords={'obsloc':np.arange(npts)})
            
        # Observation lat/lon
        if (date==str(sdate)):
            aermsk=aermsk_tmpds 
        else:
            aermsk=xa.concat((aermsk,aermsk_tmpds),dim='obsloc')
    
    total_obscounts=aermsk.obsloc.size
    aermsk=aermsk.assign_coords(obsloc=np.arange(total_obscounts))
   
    df=aermsk.to_dataframe()
 
    latbin=np.arange(-90,90+0.5*degres,degres)
    latgrd=np.arange(-90+0.5*degres,90,degres)
    lonbin=np.arange(-180,180+0.5*degres,degres)
    longrd=np.arange(-180+0.5*degres,180,degres)
    
    df['lat']=pd.cut(df['rlat'],bins=latbin,labels=latgrd)
    df['lon']=pd.cut(df['rlon'],bins=lonbin,labels=longrd)
    
    grp = df.groupby(['lat','lon']).agg({'aer_msk':['count']})
    grp=grp.rename(columns={'aer_msk':'total'})

    df_filter=((df['aer_msk']==True))
    filtered_df=df.loc[df_filter,:]   
    filtered_df=filtered_df.reset_index()
    
    filtered_df['lat']=pd.cut(filtered_df['rlat'],bins=latbin,labels=latgrd)
    filtered_df['lon']=pd.cut(filtered_df['rlon'],bins=lonbin,labels=longrd)
    filtered_grp = filtered_df.groupby(['lat','lon']).agg({'aer_msk':['count']})
    filtered_grp = filtered_grp.rename(columns={'aer_msk':'hazy'})

    outgrp=pd.concat((grp,filtered_grp),axis=1)
    
    grdds=outgrp.to_xarray()
    for var in ['total','hazy']:
        newname='%s_count'%(var)
        grdds=grdds.rename({(var,'count'):(newname)})

savedir=outputpath
if ( not os.path.exists(savedir) ):
    os.makedirs(savedir)
filename='%s/aer_msk.%.1fx%.1f.%s_%s.nc4' %(savedir,degres,degres,sdate,edate)
if (gen_file):
    grdds.to_netcdf(filename)

if (gen_plot):
    print('Generate plot')
    imgsavpath=rootpath+'/AlbanyWork/Prospectus/Experiments/HazyDA/Images/DiagFiles/gridded/2dmap/aermsk/'+area
    if ( not os.path.exists(imgsavpath) ):
        os.makedirs(imgsavpath)

    cnlvs=np.array((0.,1., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100.))
    clridx=np.array((2,3,5,7,8,9,10,12,14,16,17,18))
    clrmap=setup_cmap('precip3_16lev',clridx)
    aer_norm = mpcrs.BoundaryNorm(cnlvs,len(clridx),extend='max')
    cblb='Fraction of aerosol-affected data [%]'

    if (not gen_data):
       grdds=xa.open_dataset(filename)

    pltdata=grdds['hazy_count']/grdds['total_count']*100
    #pltdata=ds.aero_frac
    fig,ax,gl=setupax_2dmap(cornll,area,proj,lbsize=lbsize)
    set_size(axe_w,axe_h,b=0.13,l=0.05,r=0.95,t=0.95)
    cn=ax.contourf(pltdata.lon,pltdata.lat,pltdata,levels=cnlvs,cmap=clrmap,norm=aer_norm,extend='max')
    plt.colorbar(cn,ax=ax,orientation='horizontal',
                 fraction=0.045,aspect=40,pad=0.08,label=cblb)
    outname='%s/aer_msk_%.1fx%.1f.%s_%s.%s' %(imgsavpath,degres,degres,sdate,edate,ffmt)
    print(outname)
    fig.savefig(outname,dpi=quality)
    plt.close()
