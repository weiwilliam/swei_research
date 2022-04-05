#!/usr/bin/env python3
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
    rootpath='/glade/work/swei/output/images/Dataset'
    rootarch='/scratch2/BMC/gsd-fv3-dev/Shih-wei.Wei/ResearchData'
    rootgit='/glade/u/home/swei/research'
sys.path.append(rootgit+'/pyscripts/functions')
import setuparea as setarea
from plot_utils import set_size
from utils import ndate
from datetime import datetime, timedelta
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import xarray as xa
import pandas as pd
from scipy import stats

tlsize=10 ; lbsize=8
mpl.rc('axes', titlesize=tlsize,labelsize=lbsize)
mpl.rc('xtick',labelsize=lbsize)
mpl.rc('ytick',labelsize=lbsize)
mpl.rc('legend',fontsize='small')
fsave=1 ; ffmt='png' ; ptsize=4
axe_w=3 ; axe_h=3 ; quality=300


###################################################
#### UI: fhrs:24,48,72,66,60,54; yvar:tsfc,swdn ####
fhrs=24
pltvar='tsfc'
###################################################

exp0path='/glade/scratch/swei/UFS_SMOKE'
exp1path='/glade/scratch/clu/UFS_SMOKE_RR'
outputpath=rootpath+'/UFS/xyscatter'
if ( not os.path.exists(outputpath) ):
    os.makedirs(outputpath)

expnlist=['clim','rr06']
explist=['ufs.mfu0','ufs.mfu6']

xvar='aod550'
xunits=''
if (pltvar=='tsfc'):
  yvar='tmpsfc'
  yunits='K'
  ylim=[265.,335.]
#if (pltvar=='swdn'):
#  yvar='dswrf'
#  yunits='W/m2'
#  ylim=[0.,750.]

xlb='%s [%s]; %s'%(pltvar,yunits,expnlist[0])
ylb='%s [%s]; %s'%(pltvar,yunits,expnlist[1])
plot_stratified_x=1
level_x=[0.1,0.6,1.5]
cldcv_thres=0.4

sdate=2020082200
edate=2020092100
hint=24

clrlst=['tab:green','tab:blue','tab:red']
alplst=[0.6, 0.4, 0.2]

area='CONUS'# Glb, NPO, NML, TRO, SML, SPO, EAsia, NAfr
minlon, maxlon, minlat, maxlat, crosszero, cyclic=setarea.setarea(area)
print(minlat,maxlat,minlon,maxlon,crosszero, cyclic)

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
    cdate=ndate(hint,cdate)

date_counts=0
for date in dlist:
    fcstfile='phyf%.3i.nc'%(fhrs)
    infile0=exp0path+'/'+explist[0]+'.'+date+'/'+fcstfile
    infile1=exp1path+'/'+explist[1]+'.'+date+'/'+fcstfile
    print('infile0: %s' %infile0)
    print('infile1: %s' %infile1)
    if (os.path.exists(infile0) and os.path.exists(infile1)):
        ds0=xa.open_dataset(infile0)
        ds0=ds0.sel(time=ds0.time[0])
        ds1=xa.open_dataset(infile1)
        ds1=ds1.sel(time=ds1.time[0])
    else:
        print('%s is not existing'%(fcstfile))
##      d=d+1
        continue

    if (area=='Glb'):
        tmpvar0_x=ds0[xvar]
        tmpvar1_x=ds1[xvar]
        tmpvar0_y=ds0[yvar]
        tmpvar1_y=ds1[yvar]
        tmpcld0=ds0['tcdc_aveclm']
        tmpcld1=ds1['tcdc_aveclm']
    else:
        tmpvar0_x=ds0[xvar].sel(grid_xt=slice(minlon,maxlon),grid_yt=slice(maxlat,minlat))
        tmpvar1_x=ds1[xvar].sel(grid_xt=slice(minlon,maxlon),grid_yt=slice(maxlat,minlat))
        tmpvar0_y=ds0[yvar].sel(grid_xt=slice(minlon,maxlon),grid_yt=slice(maxlat,minlat))
        tmpvar1_y=ds1[yvar].sel(grid_xt=slice(minlon,maxlon),grid_yt=slice(maxlat,minlat))
        tmpcldcv0=ds0['tcdc_aveclm'].sel(grid_xt=slice(minlon,maxlon),grid_yt=slice(maxlat,minlat))
        tmpcldcv1=ds1['tcdc_aveclm'].sel(grid_xt=slice(minlon,maxlon),grid_yt=slice(maxlat,minlat))
    
    tmpvar0_x=tmpvar0_x.stack(grids=('grid_yt','grid_xt'))
    tmpvar1_x=tmpvar1_x.stack(grids=('grid_yt','grid_xt'))
    tmpvar0_y=tmpvar0_y.stack(grids=('grid_yt','grid_xt'))
    tmpvar1_y=tmpvar1_y.stack(grids=('grid_yt','grid_xt'))
    tmpcldcv0=tmpcldcv0.stack(grids=('grid_yt','grid_xt'))
    tmpcldcv1=tmpcldcv1.stack(grids=('grid_yt','grid_xt'))

    if (date_counts==0):
       var0_x=tmpvar0_x
       var1_x=tmpvar1_x
       var0_y=tmpvar0_y
       var1_y=tmpvar1_y
       cldcv0=tmpcldcv0
       cldcv1=tmpcldcv1
    else:
       var0_x=xa.concat((var0_x,tmpvar0_x),dim='grids')
       var1_x=xa.concat((var1_x,tmpvar1_x),dim='grids')
       var0_y=xa.concat((var0_y,tmpvar0_y),dim='grids')
       var1_y=xa.concat((var1_y,tmpvar1_y),dim='grids')
       cldcv0=xa.concat((cldcv0,tmpcldcv0),dim='grids')
       cldcv1=xa.concat((cldcv1,tmpcldcv1),dim='grids')

    date_counts+=1
##x: aod; y: tsfc
##0: clim; 1: rr06       
var_x=var0_y
var_y=var1_y
var_z=var1_x-var0_x

title='%s vs %s (FH%s)' %(expnlist[0],expnlist[1],fhrs)
fname='%s/%s_%s_%s_%s_fh%s.%s'%(outputpath,area,pltvar,sdate,edate,fhrs,ffmt)
print(fname) 

fig,ax=plt.subplots()
set_size(axe_w,axe_h,l=0.15,b=0.15)
if (plot_stratified_x):
   i=0
   lglst=[]
   for i in np.arange(len(level_x)):
       if (i!=len(level_x)-1):
          tmpmsk=(var_z>level_x[i])&(var_z<=level_x[i+1])
          lglst.append('aod diff (%.1f-%.1f)' %(level_x[i],level_x[i+1]))
       else:
          tmpmsk=(var1_x>level_x[i])
          lglst.append('aod diff (%.1f- )' %(level_x[i]))
       tmpmsk=(tmpmsk)&(cldcv1<=cldcv_thres)

       ax.set_xlim(xmin=ylim[0],xmax=ylim[1])
       ax.set_ylim(ymin=ylim[0],ymax=ylim[1])
       sc=ax.scatter(var_x[tmpmsk],var_y[tmpmsk],c=clrlst[i],
                     s=ptsize,alpha=alplst[i],edgecolors='none')
       stats_result = stats.linregress(var_x[tmpmsk],var_y[tmpmsk])
       print( '%s:' %(lglst[i]) )
       print( 'slope= %.2f' %(stats_result.slope) )
       print( 'intercept= %.2f' %(stats_result.intercept) )
       print( 'r-value= %.2f' %(stats_result.rvalue) )
       print( 'p-value= %.2f' %(stats_result.pvalue) )
       print( 'std err= %.2f' %(stats_result.stderr) )
       print( 'intercept std err= %.2f' %(stats_result.intercept_stderr) )
else:
   sc=ax.scatter(var_x,var_y,c='tab:blue',s=ptsize,edgecolors='none')
   lglst=['all grids']

ax.legend(lglst)
ax.set_xlabel(xlb)
ax.set_ylabel(ylb)
ax.set_title(title,loc='left')
 
if (fsave):
   fig.savefig(fname,dpi=quality)
   plt.close()
    

