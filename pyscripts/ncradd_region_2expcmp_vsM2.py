#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 21:59:35 2019

@author: weiwilliam
"""
import os,sys,platform
os_name=platform.system()
if (os_name=='Darwin'):
    rootpath='/Users/weiwilliam'
    rootarch='/Volumes/WD2TB/ResearchData'
elif (os_name=='Windows'):
    rootpath='F:\GoogleDrive_NCU\Albany'
    rootarch='F:\ResearchData'
    rootgit='F:\GitHub\swei_research'
sys.path.append(rootpath+'/AlbanyWork/Utility/Python3/functions')
import setuparea as setarea
from utils import ndate,setup_cmap
from plot_utils import setupax_2dmap,set_size
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mpcrs
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from datetime import datetime
from datetime import timedelta
#import matplotlib.dates as mdates
#from matplotlib.dates import (DAILY, DateFormatter,
#                              rrulewrapper, RRuleLocator)
#from matplotlib.ticker import MaxNLocator
#from matplotlib.gridspec import GridSpec
import time
import xarray as xa
import pandas as pd

mpl.rc('axes',titlesize=12,labelsize=12)
mpl.rc('xtick',labelsize=12)
mpl.rc('ytick',labelsize=12)
mpl.rc('legend',fontsize='x-large')
minussign=u'\u2212'
fsave=1 ; ffmt='png'
axe_w=4 ; axe_h=4 ; gllbsize=12

proj=ccrs.PlateCarree(globe=None)
grav=9.80665e+0
outpath=rootpath+'/AlbanyWork/Prospectus/Experiments/SingleCycleGSI/Images/ScatterXY'
if ( not os.path.exists(outpath) ):
   os.makedirs(outpath)

path=rootarch+'/Prospectus/SingleCycleGSI/ncdiag'
m2path=rootarch+'/common/MERRA2'
m2tag='inst3_3d_aer_Nv'
date=2020062212
plotimg='OMFdiff' # OMFdiff simBTdiff obsBTall
pltqc=0
usebc=0
ptsize=6 # Glb:2 other:10
add_insert=0
aersp='dust' # total, dust, salt, carbon, sulf

sensor='iasi_metop-a'

expname1='control'
expn1='CTL'
expname2='aerorad'
expn2='AER'

chkwvl=10.39

area='r2o1'
minlon, maxlon, minlat, maxlat, crosszero, cyclic=setarea.setarea(area)
print(area,minlat,maxlat,minlon,maxlon,crosszero,cyclic)
corll=[minlat,maxlat,minlon,maxlon]

tbclridx=[2,11,20,29,38,47,56,65,74,83,92,101,110,119,128]
tb_cmap=setup_cmap('MPL_jet',tbclridx)
tbmin=285.; tbmax=313.
tbcbticks=np.linspace(tbmin,tbmax,len(tbclridx))
tbnorm=mpcrs.BoundaryNorm(tbcbticks,len(tbclridx)+1,extend='both')

diffcidx=[2,4,6,8,10,12,14,16,18,20]
df_cmap=setup_cmap('temp_diff_18lev',diffcidx)
dfvmax=5
dfcbticks=np.linspace(-dfvmax,dfvmax,len(diffcidx)-1)
dfnorm=mpcrs.BoundaryNorm(dfcbticks,len(diffcidx),extend='both')

conc_clridx=[2,11,20,29,38,47,56,65,74,83,92,101,110,119,128]
conc_cmap=setup_cmap('MPL_jet',conc_clridx)
conc_lvs=[1e-9, 5e-9, 1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 
          5e-6, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 1e-2]
conc_norm = mpcrs.BoundaryNorm(conc_lvs,len(conc_clridx)+1,extend='both')

loop='ges' #ges,anl
if loop=='anl':
    tlstr='OMA'
elif loop=='ges':
    tlstr='OMF'

if (pltqc==-1):
    qcstr='noqc'
else:
    qcstr='qc%s'%(pltqc)

if (usebc):
   bcstr='bc'
else:
   bcstr='nobc'

raddfile='diag_'+sensor+'_'+loop+'.'+str(date)+'.nc4'
infile1=path+'/'+expname1+'/'+str(date)+'/'+raddfile
infile2=path+'/'+expname2+'/'+str(date)+'/'+raddfile

if (os.path.exists(infile1) and os.path.exists(infile2)):
    print('Processing Radfile: %s' %(raddfile))
    ds1=xa.open_dataset(infile1)
    ds2=xa.open_dataset(infile2)
    npts=int(ds1.nobs.size/ds1.nchans.size)
    nchs=ds1.nchans.size
    wavelength=1e+04/ds1.wavenumber
    usedchidx=np.where(ds1.use_flag==1)[0]
    unusedchidx=np.where(ds1.use_flag==-1)[0]
    wvldiff=abs(np.subtract(wavelength[usedchidx],chkwvl))
    chkwvlidx=usedchidx[np.where(wvldiff==wvldiff.min())[0][0]]
    print('Check wavelength: %.2f' %(wavelength[chkwvlidx]))
else:
    print('%s is not existing'%(raddfile))
    sys.exit()
        
rlat1=np.reshape(ds1.Latitude.values,(npts,nchs))[:,chkwvlidx]
rlon1=np.reshape(ds1.Longitude.values,(npts,nchs))[:,chkwvlidx]
qcflags1=np.reshape(ds1.QC_Flag.values,(npts,nchs))[:,chkwvlidx]
qcflags2=np.reshape(ds2.QC_Flag.values,(npts,nchs))[:,chkwvlidx]

locidx=np.arange(npts)

mask=~np.isnan(locidx)

if (crosszero):
    rlon1[rlon1>=maxlon]=rlon1[rlon1>=maxlon]-360.

if (area!='Glb'):
    mask=(mask)&((rlon1<maxlon)&(rlon1>minlon)&(rlat1>minlat)&(rlat1<maxlat))

if (pltqc!=-1):
    mask=(mask)&(qcflags2==pltqc)

if (plotimg=='obsBT'):
   obs1=np.reshape(ds1.Observation.values,(npts,nchs))[:,chkwvlidx]
   pltdata=obs1   
   pltidx=locidx[mask==1]
   pltnum=np.count_nonzero(mask)
   colormap=tb_cmap 
   cblbfmt='%.1f'
   clbstr='BT [K]'
   norm=tbnorm
   titlestr='Observation BT (%.2f K)' %(pltnum,np.nanmean(pltdata[pltidx]))
   fname='%s_%s_%s_ObsBT_%.2f_%s.%s'%(sensor,area,expn1,wavelength[chkwvlidx],date,ffmt)

if (plotimg=='simBTdiff'):
   simtb1=np.reshape(ds1.Simulated_Tb.values,(npts,nchs))[:,chkwvlidx]
   simtb2=np.reshape(ds2.Simulated_Tb.values,(npts,nchs))[:,chkwvlidx]
   pltdata=simtb2-simtb1
   pltidx=locidx[mask==1]
   pltnum=np.count_nonzero(mask)
   colormap=df_cmap 
   cblbfmt='%.2f'
   clbstr='BTD [K]'
   norm=dfnorm
   titlestr='%s%s%s SimBT (Mean: %.2f K)' %(expn2,minussign,expn1,np.nanmean(pltdata[pltidx]))
   fname='%s_%s_%s-%s_%s_SimBT_%.2f_%s.%s'%(sensor,area,expn2,expn1,loop,wavelength[chkwvlidx],date,ffmt)

if (plotimg=='OMFdiff'):
   if (usebc):
      omf1=np.reshape(ds1.Obs_Minus_Forecast_adjusted.values,(npts,nchs))[:,chkwvlidx]
      omf2=np.reshape(ds2.Obs_Minus_Forecast_adjusted.values,(npts,nchs))[:,chkwvlidx]
   else:
      omf1=np.reshape(ds1.Obs_Minus_Forecast_unadjusted.values,(npts,nchs))[:,chkwvlidx]
      omf2=np.reshape(ds2.Obs_Minus_Forecast_unadjusted.values,(npts,nchs))[:,chkwvlidx]
   pltdata=omf2-omf1
   pltidx=locidx[mask==1]
   rmsd1=np.sqrt(np.nanmean(np.square(omf1[pltidx])))
   rmsd2=np.sqrt(np.nanmean(np.square(omf2[pltidx])))
   rmsd_df=rmsd2-rmsd1
   print("%s RMSD1: %.3f" %(area,rmsd1))
   print("%s RMSD2: %.3f" %(area,rmsd2))
   print("%s RMSD(2-1): %.3f" %(area,rmsd_df))
   titlestr='%s%s%s' %(expn2,minussign,expn1)
   xlbstr='OMF Difference [K]'
   ylbstr='%s column mass density [$\mathrm{kg\cdot m^{-2}}$]' %(aersp.capitalize())
   # titlestr='%s%s%s %s (Mean: %.2f K RMS Diff: %.2f K)' %(expn2,minussign,expn1,tlstr,np.nanmean(pltdata[pltidx]),rmsd_df)
   fname='%s_%s_%s-%s_%s_%.2f_%s_%s.%s.%s'%(sensor,area,expn2,expn1,tlstr,wavelength[chkwvlidx],bcstr,qcstr,date,ffmt)

# MERRA2
if (aersp=='dust'):
    varlst=['DU001','DU002','DU003','DU004','DU005']
elif (aersp=='salt'):
    varlst=['SS001','SS002','SS003','SS004','SS005']
elif (aersp=='carbon'):
    varlst=['OCPHILIC','OCPHOBIC','BCPHILIC','BCPHOBIC']
elif (aersp=='sulf'):
    varlst=['SO4']
elif (aersp=='total'):
    varlst=['DU001','DU002','DU003','DU004','DU005',
            'SS001','SS002','SS003','SS004','SS005',
            'OCPHILIC','OCPHOBIC','BCPHILIC','BCPHOBIC',
            'SO4']
cnlvs=[4e-4, 8e-4, 1.6e-3]
yy=str(date)[:4] ; mm=str(date)[4:6]
dd=str(date)[6:8]; hh=str(date)[8:10]
pdy=str(date)[:8]
cdate=datetime(int(yy),int(mm),int(dd),int(hh))
m2file=m2path+'/'+yy+'/'+mm+'/MERRA2_400.'+m2tag+'.'+pdy+'.nc4'
if (os.path.exists(m2file)):
    m2=xa.open_dataset(m2file)
    m2=m2.sel(time=cdate)
    if (area!='Glb'):
       m2=m2.sel(lon=slice(minlon-5,maxlon+5),lat=slice(minlat-5,maxlat+5))
    delp=m2.DELP
    kgkg_kgm2=delp/grav

    for var in varlst:
        tmp=m2[var]*kgkg_kgm2
        if (var==varlst[0]):
            conc=tmp
        else:
            conc=xa.concat((conc,tmp),dim='bins')
    cmass=conc.sum(dim=('bins','lev'))

for i in np.arange(pltidx.size):
    intidx=pltidx[i]
    tmp=cmass.interp(lat=rlat1[intidx],lon=rlon1[intidx],method='cubic')
    if (i==0):
        cmass_at_obs=tmp
    else:
        cmass_at_obs=xa.concat((cmass_at_obs,tmp),dim='obsloc')

fig,ax=plt.subplots()
set_size(axe_w,axe_h,l=0.15,b=0.15)
sc=ax.scatter(pltdata[pltidx],cmass_at_obs,c='tab:blue',s=ptsize)
ax.set_xlabel(xlbstr)
ax.set_ylabel(ylbstr)
ax.set_title(titlestr)
ax.spines.right.set_visible(False)
ax.spines.top.set_visible(False)
ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

if (add_insert):
    log_ax = fig.add_axes([0.6, 0.15, 0.3, 0.3])
    log_ax.set_xlim(-0.1,2)
    log_ax.set_ylim(1e-6,1e-4)
    log_ax.set_yscale('log')
    logsc=log_ax.scatter(pltdata[pltidx],cmass_at_obs,c='tab:blue',s=ptsize)
    log_ax.xaxis.tick_top()
    log_ax.yaxis.tick_right()
    log_ax.xaxis.set_tick_params(rotation=-45)

if (fsave):
   fig.savefig(outpath+'/'+fname,dpi=200)
   plt.close()

