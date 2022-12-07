#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 18 12:41:00 2018

@author: weiwilliam
"""
import os,sys,platform
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
import setuparea as setarea
from plot_utils import set_size
from utils import gen_eqs_by_stats,ndate
from time import time
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.dates import (DAILY, DateFormatter,
                              rrulewrapper, RRuleLocator)
import matplotlib.cbook as cbook
import numpy as np
import scipy.io as io
from scipy import stats
from datetime import datetime, timedelta
import xarray as xa
import seaborn as sns
import pandas as pd
from gsi_ncdiag import read_cnv_ncdiag

axe_w=6 ; axe_h=6
mpl.rc('axes',titlesize=16,labelsize=16)
mpl.rc('xtick',labelsize=16)
mpl.rc('ytick',labelsize=16)
mpl.rc('legend',fontsize='x-large')
fsave=1; ffmt='png'; quality=300

sfctype_list=['180','181','182','183','187','189','196','197','198','199']
varlist=['sst'] #['ps','sst','gps','q','t','uv','tcp']
unitlist=['K'] #['mb','K','%','g/kg','K','m/s','mb']
colorlist=['b','r']
markerlist=['*','o']

explist=['hazyda_ctrl','hazyda_aero']
expnlist=['CTL','AER']

sdate=2020062100
edate=2020062518
hint=6
bufrtype='all'
loop='ges' #ges,anl
useqc=1
usebc=0
area='r2o7'# Glb, NPO, NML, TRO, SML, SPO, EAsia, NAfr

inputpath=rootarch
outputpath=rootpath+'/AlbanyWork/Prospectus/Experiments/HazyDA/Images/DiagFiles/conv/scatter'
imgsavpath=outputpath+'/'+area
if ( not os.path.exists(imgsavpath) ):
   os.makedirs(imgsavpath)

if (loop=='anl'):
    tlstr='OMA'
elif (loop=='ges'):
    tlstr='OMF'
if (useqc):
    qcflg='qc'
else:
    qcflg='noqc'

minlon, maxlon, minlat, maxlat, crosszero, cyclic=setarea.setarea(area)
print(minlat,maxlat,minlon,maxlon,crosszero, cyclic)
if (cyclic):
    minlon=-180.
    maxlon=180.
else:
    minlon=(minlon+180)%360-180
    maxlon=(maxlon+180)%360-180

def scatterplot_x2y2(xdata1,xdata2,xlb,ydata1,ydata2,ylb,stat1,stat2,lglst,clrlst,mkrlst,title,fname,fig=None,ax=None):
    from utils import gen_eqs_by_stats
    if not fig: fig=plt.gcf()
    if not  ax: ax=plt.gca()
    ax.scatter(xdata1,ydata1,c=clrlst[0],marker=mkrlst[0],alpha=0.5,edgecolors='None')
    ax.scatter(xdata2,ydata2,c=clrlst[1],marker=mkrlst[1],alpha=0.5,edgecolors='None')

    vmin=np.floor(np.array((xdata1.min(),xdata2.min(),ydata1.min(),ydata2.min())).min())
    vmax= np.ceil(np.array((xdata1.max(),xdata2.max(),ydata1.max(),ydata2.max())).max())
    ax.set_xlim(vmin,vmax)
    ax.set_ylim(vmin,vmax)

    regeqs1=gen_eqs_by_stats(stat1)
    regeqs2=gen_eqs_by_stats(stat2)
    lglst[0]='%s (%s) $R^2=$ %.2f' %(lglst[0],regeqs1,np.square(stat1.rvalue))
    lglst[1]='%s (%s) $R^2=$ %.2f' %(lglst[1],regeqs2,np.square(stat2.rvalue))
    ax.legend(lglst)

    ax.plot(np.arange(vmin,vmax+.5),np.arange(vmin,vmax+.5),color='grey',linewidth=0.8)
    ax.plot(xdata1,stat1.slope*xdata1+stat1.intercept,color=clrlst[0])
    ax.plot(xdata2,stat2.slope*xdata2+stat2.intercept,color=clrlst[1])
    if (np.arange(vmin,vmax+.5).size>10 or
        np.arange(vmin,vmax+.5).size<4):
       xyticknum=8
    else:
       xyticknum=np.arange(vmin,vmax+.5).size
    ax.locator_params(axis='x',nbins=xyticknum)
    ax.locator_params(axis='y',nbins=xyticknum)
    ax.set_xlabel(xlb)
    ax.set_ylabel(ylb)
    ax.set_title(title,loc='left')
    fig.savefig(fname,dpi=300)
    plt.close()
    return
    #sns.regplot(x=xdata1,y=ydata1,color=clrlst[0],marker=mkrlst[0],scatter_kws={"alpha":0.7},ax=ax)
    #sns.regplot(x=xdata2,y=ydata2,color=clrlst[1],marker=mkrlst[1],scatter_kws={"alpha":0.7},ax=ax)

date1 = pd.to_datetime(sdate,format='%Y%m%d%H')
date2 = pd.to_datetime(edate,format='%Y%m%d%H')
delta = timedelta(hours=6)
dates = pd.date_range(start=date1, end=date2, freq=delta)

rule = rrulewrapper(DAILY, byhour=6, interval=5)
loc = RRuleLocator(rule)
formatter = DateFormatter('%Y%h %n %d %Hz')

# Calculate how many cases are going to do statistic.
tnum=0
dlist=[]
cdate=sdate
while (cdate<=edate):
    dlist.append(str(cdate))
    tnum=tnum+1
    cdate=ndate(hint,cdate)
    
print('Total cases number is %d' % tnum )
ptop=np.array((1000.,900.,800.,600.,400.,300.,250.,200.,150.,100.,50.,0.))
pbot=np.array((1200.,1000.,900.,800.,600.,400.,300.,250.,200.,150.,100.,50.))
znum=ptop.size

uidx=0
for var in varlist:
    unit=unitlist[uidx]
    sfcflag=(var=='ps' or var=='sst' or var=='tcp'
             or bufrtype in sfctype_list)
    if (sfcflag):
        d=0
        for date in dlist:
            cnvdfile='diag_conv_'+var+'_'+loop+'.'+date+'.nc4'
            print('File: %s' % cnvdfile)
            infile1=inputpath+'/'+explist[0]+'/'+date+'/'+cnvdfile
            infile2=inputpath+'/'+explist[1]+'/'+date+'/'+cnvdfile
            if (os.path.exists(infile1) and os.path.exists(infile2)):
                print('Processing Cnvfile: %s' %(cnvdfile))
                #ds1=xa.open_dataset(infile1)
                #ds2=xa.open_dataset(infile2)
            else:
                print('%s is not existing'%(cnvdfile))
                d=d+1
                continue

            tmpds1=read_cnv_ncdiag(infile1)
            tmpds2=read_cnv_ncdiag(infile2)
            
            if (date==dlist[0]):
                scatda1=tmpds1
                scatda2=tmpds2
            else:
                scatda1=xa.concat((scatda1,tmpds1),dim='obsloc')
                scatda2=xa.concat((scatda2,tmpds2),dim='obsloc')

        mask1=~np.isnan(scatda1.obsloc)
        mask2=~np.isnan(scatda2.obsloc)

        if (area!='Glb'):
            mask1=(mask1)&((scatda1.rlon<maxlon)&(scatda1.rlon>minlon)
                           &(scatda1.rlat>minlat)&(scatda1.rlat<maxlat))
            mask2=(mask2)&((scatda2.rlon<maxlon)&(scatda2.rlon>minlon)
                           &(scatda2.rlat>minlat)&(scatda2.rlat<maxlat))
        
        if (useqc):
            mask1=(mask1)&(scatda1.qcflag==1)
            mask2=(mask2)&(scatda2.qcflag==1)

        if (bufrtype!='all'):
            mask1=(mask1)&(scatda1.obstype==int(bufrtype))
            mask2=(mask2)&(scatda2.obstype==int(bufrtype))
        
        cnts1=np.count_nonzero(mask1)
        cnts2=np.count_nonzero(mask2)
        
        pltobs1=scatda1.obs[mask1]
        pltobs2=scatda2.obs[mask2]
        if (usebc):
           pltmdl1=(scatda1.obs-scatda1.omb_bc)[mask1]
           pltmdl2=(scatda2.obs-scatda2.omb_bc)[mask2]
        else:
           pltmdl1=(scatda1.obs-scatda1.omb_nbc)[mask1]
           pltmdl2=(scatda2.obs-scatda2.omb_nbc)[mask2]

        pltmean1=np.mean(pltobs1-pltmdl1)
        pltmean2=np.mean(pltobs2-pltmdl2)
        pltrms1=np.sqrt(np.mean(np.square(pltobs1-pltmdl1)))
        pltrms2=np.sqrt(np.mean(np.square(pltobs2-pltmdl2)))

        stat_res1=stats.linregress(pltobs1,pltmdl1)
        stat_res2=stats.linregress(pltobs2,pltmdl2)
        
        xlb='Observed %s [%s]' %(var.upper(),unit)
        if (loop=='ges'):
           ylb='First-guess %s [%s]'%(var.upper(),unit)
        elif (loop=='anl'):
           ylb='Analyzed %s [%s]'%(var.upper(),unit)

        title='RMS %s:%.2f %s:%.2f Mean %s:%.2f %s:%.2f' %(expnlist[0], pltrms1,expnlist[1], pltrms2,
                                                           expnlist[0],pltmean1,expnlist[1],pltmean2)
        fname='%s/%s_%s_%s_%s_bufr%s.%s_%s.%s_%s.%s'%(imgsavpath,area,var.upper(),tlstr,qcflg,
                                                   bufrtype,cnts1,cnts2,sdate,edate,ffmt)
        print(fname,flush=1)

        fig,ax=plt.subplots()
        set_size(axe_w,axe_h,ax=ax,l=0.15)
        scatterplot_x2y2(pltobs1,pltobs2,xlb,pltmdl1,pltmdl2,ylb,stat_res1,stat_res2,expnlist,
                         colorlist,markerlist,title,fname,ax=ax)
         
    else:
        for z in np.arange(1):
            d=0
            for date in dlist:
                cnvdfile='diag_conv_'+var+'_'+loop+'.'+date+'.nc4'
                print('File: %s' % cnvdfile)
                infile1=path+'/'+explist[0]+'/'+date+'/'+cnvdfile
                infile2=path+'/'+explist[1]+'/'+date+'/'+cnvdfile
                if (os.path.exists(infile1) and os.path.exists(infile2)):
                    print('Processing Cnvfile: %s' %(cnvdfile))
                    ds1=xa.open_dataset(infile1)
                    ds2=xa.open_dataset(infile2)
                    try:
                        omg_mean
                    except NameError:
                        omg_mean=np.zeros((tnum,znum,2),dtype='float')
                        omg_rmsq=np.zeros_like(omg_mean)
                        omg_mean[:,:,:]=np.nan
                        omg_rmsq[:,:,:]=np.nan
                else:
                    print('%s is not existing'%(cnvdfile))
                    d=d+1
                    continue
        
                rlat1=ds1.obsdata[2,:]
                rlon1=ds1.obsdata[3,:]
                rlat2=ds2.obsdata[2,:]
                rlon2=ds2.obsdata[3,:]
        
                if (crosszero):
                    rlon1[rlon1>=maxlon]=rlon1[rlon1>=maxlon]-360.
                    rlon2[rlon2>=maxlon]=rlon2[rlon2>=maxlon]-360.
        
                pres1=ds1.obsdata[5,:]
                iuse1=ds1.obsdata[11,:]
                pres2=ds2.obsdata[5,:]
                iuse2=ds2.obsdata[11,:]
        
                mask1=(ds1.vartype==var)
                mask2=(ds2.vartype==var)
        
                if (area!='Glb'):
                    mask1=(mask1)&((rlon1<maxlon)&(rlon1>minlon)&(rlat1>minlat)&(rlon1<maxlat))
                    mask2=(mask2)&((rlon2<maxlon)&(rlon2>minlon)&(rlat2>minlat)&(rlon2<maxlat))
        
                if (useqc):
                    mask1=(mask1)&(iuse1==1)
                    mask2=(mask2)&(iuse2==1)
            
                if (var=='uv'):
                    dpar1=ds1.u_Obs_Minus_Forecast_adjusted
                    dpar2=ds2.u_Obs_Minus_Forecast_adjusted
                else:
                    dpar1=ds1.Obs_Minus_Forecast_adjusted
                    dpar2=ds2.Obs_Minus_Forecast_adjusted
                
                zmask1=(mask1)&((pres1<pbot[z])&(pres1>ptop[z]))
                zmask2=(mask2)&((pres2<pbot[z])&(pres2>ptop[z]))
                
                zdpar1=xa.where(~zmask1,np.nan,dpar1)
                zdpar2=xa.where(~zmask2,np.nan,dpar2)
        
                try:
                    scatda1
                except NameError:
                    pltdata1=zdpar1 ; pltdata2=zdpar2
                    scatda1=pltdata1 ; lon1=rlon1 ; lat1=rlat1
                    scatda2=pltdata2 ; lon2=rlon2 ; lat2=rlat2
                else:
                    pltdata1=zdpar1 ; pltdata2=zdpar2
                    scatda1=xa.concat((scatda1,pltdata1),dim='obsloc')
                    lon1=xa.concat((lon1,rlon1),dim='obsloc')
                    lat1=xa.concat((lat1,rlat1),dim='obsloc')
                    scatda2=xa.concat((scatda2,pltdata2),dim='obsloc')
                    lon2=xa.concat((lon2,rlon2),dim='obsloc')
                    lat2=xa.concat((lat2,rlat2),dim='obsloc')
                
                vmax=np.max((np.nanmax(pltdata1),np.nanmax(pltdata2)))
            
                pltnum=np.count_nonzero(~np.isnan(pltdata1))
                pltrms=np.sqrt(np.nanmean(np.square(pltdata1)))
                title1=('%s %s %s (Mean:%.2f,RMS:%.2f,OBS#:%i)'
                                   %(expnlist[0],var.upper(),tlstr,pltdata1.mean(),pltrms,pltnum))
                fname1='%s/%s_%s_%s_P%i_%s_%s.%s.%s'%(outpath,area,expnlist[0],var.upper(),
                                                   pbot[z],tlstr,qcflg,date,ffmt)
                scatterplot(axe_w,axe_h,rlon1,rlat1,pltdata1,pltrms,vmax,title1,fname1)
            
                pltnum=np.count_nonzero(~np.isnan(pltdata2))
                pltrms=np.sqrt(np.nanmean(np.square(pltdata2)))
                title2=('%s %s %s (Mean:%.2f,RMS:%.2f,OBS#:%i)'
                                   %(expnlist[1],var.upper(),tlstr,pltdata2.mean(),pltrms,pltnum))
                fname2='%s/%s_%s_%s_P%i_%s_%s.%s.%s'%(outpath,area,expnlist[1],var.upper(),
                                                   pbot[z],tlstr,qcflg,date,ffmt)
                scatterplot(fsize,rlon2,rlat2,pltdata2,pltrms,vmax,title2,fname2)
        
                pltdata=pltdata1-pltdata2
                pltnum=np.count_nonzero(~np.isnan(pltdata))
                pltrms=(np.sqrt(np.nanmean(np.square(pltdata1)))-
                        np.sqrt(np.nanmean(np.square(pltdata2))))
                vmax=np.nanmax(pltdata)
                title3=('%s-%s %s %s (Mean:%.2f,RMS:%.2f,OBS#:%i)' 
                                      %(expnlist[0],expnlist[1],var.upper(),tlstr,
                                        pltdata.mean(),pltrms,pltnum))
                fname3='%s/%s_%s-%s_%s_P%i_%s_%s.%s.%s'%(outpath,area,expnlist[0],expnlist[1],
                                                      var.upper(),pbot[z],tlstr,qcflg,date,ffmt)
                scatterplot(axe_w,axe_h,rlon2,rlat2,pltdata,pltrms,vmax,title3,fname3)
    
            vmax=np.max((np.nanmax(scatda1),np.nanmax(scatda2)))
            
            pltnum=np.count_nonzero(~np.isnan(scatda1))
            pltrms=np.sqrt(np.nanmean(np.square(scatda1)))
            title1=('%s %s %s (Mean:%.2f,RMS:%.2f,OBS#:%i)' 
                    %(expnlist[0],var.upper(),tlstr,scatda1.mean(),pltrms,pltnum))
            fname1='%s/%s_%s_%s_P%i_%s_%s.%s'%(outpath,area,expnlist[0],var.upper(),
                                               pbot[z],tlstr,qcflg,ffmt)
            scatterplot(axe_w,axe_h,lon1,lat1,scatda1,pltrms,vmax,title1,fname1)
            
            pltnum=np.count_nonzero(~np.isnan(scatda2))
            pltrms=np.sqrt(np.nanmean(np.square(scatda2)))
            title2=('%s %s %s (Mean:%.2f,RMS:%.2f,OBS#:%i)'
                    %(expnlist[1],var.upper(),tlstr,scatda2.mean(),pltrms,pltnum))
            fname2='%s/%s_%s_%s_P%i_%s_%s.%s'%(outpath,area,expnlist[1],var.upper(),
                                               pbot[z],tlstr,qcflg,ffmt)
            scatterplot(axe_w,axe_h,lon2,lat2,scatda2,pltrms,vmax,title2,fname2)
        
            pltdata=scatda1-scatda2
            pltnum=np.count_nonzero(~np.isnan(pltdata))
            pltrms=(np.sqrt(np.nanmean(np.square(scatda1)))-
                    np.sqrt(np.nanmean(np.square(scatda2))))
            vmax=np.nanmax(pltdata)
            title3=('%s-%s %s %s (Mean:%.2f,RMS:%.2f,OBS#:%i)' 
                    %(expnlist[0],expnlist[1],var.upper(),tlstr,
                      pltdata.mean(),pltrms,pltnum))
            fname3='%s/%s_%s-%s_%s_P%i_%s_%s.%s'%(outpath,area,expnlist[0],expnlist[1],
                                                  var.upper(),pbot[z],tlstr,qcflg,ffmt)
            scatterplot(axe_w,axe_h,lon2,lat2,pltdata,pltrms,vmax,title3,fname3)

    uidx=uidx+1
            
    

    

