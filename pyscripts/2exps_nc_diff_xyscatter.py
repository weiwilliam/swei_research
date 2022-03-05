#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 18 12:41:00 2018

@author: weiwilliam
"""
import os,sys
machine='Cheyenne'
if (machine=='MBP'):
    rootpath='/Users/weiwilliam'
    rootarch='/Volumes/WD2TB/ResearchData'
elif (machine=='Desktop'):
    rootpath='F:\GoogleDrive_NCU\Albany'
    rootarch='F:\ResearchData'
    rootgit='F:\GitHub\swei_research'
elif (machine=='S4'):
    rootpath='/data/users/swei'
    rootarch='/data/users/swei/ResearchData'
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
from utils import ndate
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime, timedelta
import xarray as xa
import pandas as pd

tlsize=12 ; lbsize=10
mpl.rc('axes', titlesize=tlsize,labelsize=lbsize)
mpl.rc('xtick',labelsize=lbsize)
mpl.rc('ytick',labelsize=lbsize)
mpl.rc('legend',fontsize='small')
fsave=1 ; ffmt='png' ; ptsize=4
axe_w=3 ; axe_h=3 ; quality=300
minussign=u'\u2212'

exp0path='/glade/scratch/clu/UFS_SMOKE'
exp1path='/glade/scratch/clu/UFS_SMOKE_RR'

outputpath=rootpath+'/UFS_SMOKE/xyscatter'
if ( not os.path.exists(outputpath) ):
    os.makedirs(outputpath)

xvar='aod550'
xunits=''
xlb='AOD Diff.'
yvar='tmp2m' # tmpsfc
yunits='K'
ylb='2m Temperature Diff. [%s]'%(yunits)
plot_stratified_x=1
level_x=[0.1,0.4,0.8,1.2]

expset=1
if (expset==1):
   explist=['ufs.mfu0','ufs.mfu6']
   expnlist=['CLM','RR']
elif (expset==2):
   explist=['prctrl','praero']
   expnlist=['CTL','CAER']

sdate=2020091500
edate=2020091500
hint=24
fhrs=24
clrlst=['tab:blue','tab:green','tab:orange','tab:red']

area='NAmer'# Glb, NPO, NML, TRO, SML, SPO, EAsia, NAfr
minlon, maxlon, minlat, maxlat, crosszero, cyclic=setarea.setarea(area)
print(minlat,maxlat,minlon,maxlon,crosszero, cyclic)
#if (area=='Glb'):
#   minlon=-180. ; maxlon=180.
#else:
#   minlon=(minlon+180)%360-180
#   maxlon=(maxlon+180)%360-180

syy=int(str(sdate)[:4]); smm=int(str(sdate)[4:6])
sdd=int(str(sdate)[6:8]); shh=int(str(sdate)[8:10])
eyy=int(str(edate)[:4]); emm=int(str(edate)[4:6])
edd=int(str(edate)[6:8]); ehh=int(str(edate)[8:10])

date1 = datetime(syy,smm,sdd,shh)
date2 = datetime(eyy,emm,edd,ehh)
delta = timedelta(hours=hint)
dates = pd.date_range(start=date1, end=date2, freq=delta)

# Calculate how many cases are going to do statistic.
tnum=0
dlist=[]
cdate=sdate
while (cdate<=edate):
    dlist.append(str(cdate))
    tnum=tnum+1
    cdate=ndate(6,cdate)

date_counts=0
for date in dlist:
    fcstfile='phyf%.3i.nc'%(fhrs)
    print('File: %s' % fcstfile)
    infile0=exp0path+'/'+explist[0]+'.'+date+'/'+fcstfile
    infile1=exp1path+'/'+explist[1]+'.'+date+'/'+fcstfile
    if (os.path.exists(infile0) and os.path.exists(infile1)):
        print('Processing Fcst file: %s' %(fcstfile))
        ds0=xa.open_dataset(infile0)
        ds0=ds0.sel(time=ds0.time[0])
        ds1=xa.open_dataset(infile1)
        ds1=ds1.sel(time=ds1.time[0])
    else:
        print('%s is not existing'%(fcstfile))
        d=d+1
        continue

    if (area=='Glb'):
        tmpvar0_x=ds0[xvar]
        tmpvar1_x=ds1[xvar]
        tmpvar0_y=ds0[yvar]
        tmpvar1_y=ds1[yvar]
    else:
        tmpvar0_x=ds0[xvar].sel(grid_xt=slice(minlon,maxlon),grid_yt=slice(maxlat,minlat))
        tmpvar1_x=ds1[xvar].sel(grid_xt=slice(minlon,maxlon),grid_yt=slice(maxlat,minlat))
        tmpvar0_y=ds0[yvar].sel(grid_xt=slice(minlon,maxlon),grid_yt=slice(maxlat,minlat))
        tmpvar1_y=ds1[yvar].sel(grid_xt=slice(minlon,maxlon),grid_yt=slice(maxlat,minlat))
    
    tmpvar0_x=tmpvar0_x.stack(grids=('grid_yt','grid_xt'))
    tmpvar1_x=tmpvar1_x.stack(grids=('grid_yt','grid_xt'))
    tmpvar0_y=tmpvar0_y.stack(grids=('grid_yt','grid_xt'))
    tmpvar1_y=tmpvar1_y.stack(grids=('grid_yt','grid_xt'))

    if (date_counts==0):
       var0_x=tmpvar0_x
       var1_x=tmpvar1_x
       var0_y=tmpvar0_y
       var1_y=tmpvar1_y
    else:
       var0_x=xa.concat((var0_x,tmpvar0_x),dim='grids')
       var1_x=xa.concat((var1_x,tmpvar1_x),dim='grids')
       var0_y=xa.concat((var0_y,tmpvar0_y),dim='grids')
       var1_y=xa.concat((var1_y,tmpvar1_y),dim='grids')
       
    var_x=var0_x-var1_x
    var_y=var0_y-var1_y

    title='%s%s%s' %(expnlist[0],minussign,expnlist[1])
    fname='%s/%s_%s_%s_%s_%s.%s_%s.%s'%(outputpath,area,expnlist[0],expnlist[1],
                                        xvar,yvar,sdate,edate,ffmt)
    print(fname)
    
    fig,ax=plt.subplots()
    set_size(axe_w,axe_h,l=0.15,b=0.15)
    if (plot_stratified_x):
       i=0
       lglst=[]
       for i in np.arange(len(level_x)):
           if (i!=len(level_x)-1):
              tmpmsk=(var1_x>level_x[i])&(var1_x<=level_x[i+1])
              lglst.append('%s (%.1f-%.1f)' %(xvar,level_x[i],level_x[i+1]))
           else:
              tmpmsk=(var1_x>level_x[i])
              lglst.append('%s (%.1f- )' %(xvar,level_x[i]))
           sc=ax.scatter(var_x[tmpmsk],var_y[tmpmsk],c=clrlst[i],
                         s=ptsize,alpha=0.6,edgecolors='none')
    else:
       sc=ax.scatter(var_x,var_y,c='tab:blue',s=ptsize,edgecolors='none')
       lglst=['all grids']

    ax.legend(lglst)
    ax.set_xlabel(xlb)
    ax.set_ylabel(ylb)
    ax.set_title(title,loc='left')
    #ax.spines.right.set_visible(False)
    #ax.spines.top.set_visible(False)
    #ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
     
    if (fsave):
       fig.savefig(fname,dpi=quality)
       plt.close()
    
    date_counts+=1
