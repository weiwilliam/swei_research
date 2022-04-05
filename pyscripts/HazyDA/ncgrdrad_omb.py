#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 21:59:35 2019

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
    rootpath='/data/users/swei/AlbanyWork/Prospectus/Experiments/HazyDA/Images'
    rootarch='/data/users/swei'
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
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mpcrs
import xarray as xa
import pandas as pd

import setuparea as setarea
from plot_utils import setupax_2dmap, plt_x2y, set_size
from utils import setup_cmap
import cartopy.crs as ccrs

tlsize=12 ; txsize=10
mpl.rc('axes', titlesize=tlsize,labelsize=txsize)
mpl.rc('xtick',labelsize=txsize)
mpl.rc('ytick',labelsize=txsize)
mpl.rc('legend',fontsize='small')
fsave=1 ; ffmt='png' ; ptsize=4
axe_w=6 ; axe_h=3 ; quality=300
tkfreq=2
minussign=u'\u2212'

# Projection setting
proj=ccrs.PlateCarree(globe=None)

sensorlist=['airs_aqua','amsua_aqua','amsua_metop-a','amsua_n15','amsua_n18',
            'amsua_n19','atms_npp','avhrr_metop-a','avhrr_n18','cris_npp','gmi_gpm',
            'hirs4_metop-a','hirs4_metop-b','hirs4_n19','iasi_metop-a','iasi_metop-b',
            'mhs_metop-a','mhs_metop-b','mhs_n18','mhs_n19','saphir_meghat',
            'seviri_m08','seviri_m10','sndrd1_g15','sndrd2_g15','sndrd3_g15',
            'sndrd4_g15','ssmis_f17','ssmis_f18']
#hsensorlist=['airs_aqua','iasi_metop-a','iasi_metop-b','cris_npp']
#lsensor1list=['hirs4_metop-a','hirs4_metop-b','hirs4_n19']
#lsensor2list=['sndrd1_g15','sndrd2_g15','sndrd3_g15','sndrd4_g15']
#lsensor3list=['avhrr_metop-a','avhrr_n18','seviri_m08','seviri_m10']

pltvar='omb_bc_mean'
varlb='OMB w/ BC'
sdate=2020061000
edate=2020071018
hint=6

sensor='iasi_metop-a'
chkwvn=962.5
degres=1
units='K'

explist=np.array(['hazyda_ctrl','hazyda_aero'])
expnlist=['CTL','AER']

area='Glb'
minlon, maxlon, minlat, maxlat, crosszero, cyclic=setarea.setarea(area)
print(minlat,maxlat,minlon,maxlon,crosszero,cyclic)
if (area=='Glb'):
   minlon=-180. ; maxlon=180.
else:
   minlon=(minlon+180)%360-180
   maxlon=(maxlon+180)%360-180
cornerll=[minlat,maxlat,minlon,maxlon]

loop='ges' #ges,anl
if loop=='anl':
    tlstr='OMA'
elif loop=='ges':
    tlstr='OMF'

useqc=1
if (useqc):
   qcflg='qc'
else:
   qcflg='noqc'

cbori='horizontal' #vertical, horizontal
if (cbori=='vertical'):
   cb_frac=0.025
   cb_pad=0.06
elif (cbori=='horizontal'):
   cb_frac=0.04
   cb_pad=0.1

    
inpath=rootarch+'/archive/HazyDA/gridded_diag'
outputpath=rootpath+'/DiagFiles/gridded'
imgsavpath=outputpath+'/2dmap/'+area
if ( not os.path.exists(imgsavpath) ):
   os.makedirs(imgsavpath)

grdfile0='%s/%s_%s_%s_%s_omb_%.1fx%.1f.%s_%s.nc' %(inpath,expnlist[0],sensor,loop,qcflg,degres,degres,sdate,edate)
grdfile1='%s/%s_%s_%s_%s_omb_%.1fx%.1f.%s_%s.nc' %(inpath,expnlist[1],sensor,loop,qcflg,degres,degres,sdate,edate)

ds0=xa.open_dataset(grdfile0)
ds1=xa.open_dataset(grdfile1)

pltda0=ds0[pltvar].sel(wavenumber=chkwvn)
pltda1=ds1[pltvar].sel(wavenumber=chkwvn)

lv_max=np.ceil( np.max( (np.nanquantile(abs(pltda0),0.95),np.nanquantile(abs(pltda1),0.95)) ))
lv_int=lv_max/10
tcks=np.arange(-lv_max,lv_max+lv_int,lv_int)
lvs=tcks
clridx=[]
for idx in np.linspace(2,254,tcks.size):
    clridx.append(int(idx))
clrmap=setup_cmap('BlueYellowRed',clridx)
norm = mpcrs.BoundaryNorm(lvs,len(clridx)+1,extend='both')

difflv_max=np.ceil( np.nanquantile(abs(pltda1-pltda0),0.95) )
difflv_int=difflv_max/20
difftcks=np.arange(0,difflv_max+difflv_int,difflv_int)
diff_lvs=difftcks
clridx=[]
for idx in np.linspace(2,254,difftcks.size):
    clridx.append(int(idx))
diffclrmap=setup_cmap('WhiteYellowOrangeRed',clridx)
diffnorm = mpcrs.BoundaryNorm(diff_lvs,len(clridx)+1,extend='both')

fig,ax,gl=setupax_2dmap(cornerll,area,proj,lbsize=txsize)
set_size(axe_w,axe_h,b=0.15,l=0.05,r=0.95)
pltdata=pltda0
cblabel='%s %s [%s]' %(expnlist[0],varlb,units)
cn=ax.contourf(pltdata.lon,pltdata.lat,pltdata,levels=lvs,
               cmap=clrmap,norm=norm,extend='both')
plt.colorbar(cn,orientation=cbori,fraction=cb_frac,
             pad=cb_pad,ticks=lvs[::tkfreq],aspect=40,label=cblabel)
outname='%s/%s_%s.%s_%s.%s' %(imgsavpath,expnlist[0],pltvar,sdate,edate,ffmt)
if (fsave):
   print(outname)
   fig.savefig(outname,dpi=quality)
   plt.close()

fig,ax,gl=setupax_2dmap(cornerll,area,proj,lbsize=txsize)
set_size(axe_w,axe_h,b=0.15,l=0.05,r=0.95)
pltdata=pltda1
cblabel='%s %s [%s]' %(expnlist[1],varlb,units)
cn=ax.contourf(pltdata.lon,pltdata.lat,pltdata,levels=lvs,
               cmap=clrmap,norm=norm,extend='both')
plt.colorbar(cn,orientation=cbori,fraction=cb_frac,
             pad=cb_pad,ticks=lvs[::tkfreq],aspect=40,label=cblabel)
outname='%s/%s_%s.%s_%s.%s' %(imgsavpath,expnlist[1],pltvar,sdate,edate,ffmt)
if (fsave):
   print(outname)
   fig.savefig(outname,dpi=quality)
   plt.close()

fig,ax,gl=setupax_2dmap(cornerll,area,proj,lbsize=txsize)
set_size(axe_w,axe_h,b=0.15,l=0.05,r=0.95)
pltdata=pltda1-pltda0
cblabel='%s%s%s %s [%s]' %(expnlist[1],minussign,expnlist[0],varlb,units)
cn=ax.contourf(pltdata.lon,pltdata.lat,pltdata,levels=diff_lvs,
               cmap=diffclrmap,norm=diffnorm,extend='both')
plt.colorbar(cn,orientation=cbori,fraction=cb_frac,
             pad=cb_pad,ticks=diff_lvs[::tkfreq],aspect=40,label=cblabel)
outname='%s/%s-%s_%s.%s_%s.%s' %(imgsavpath,expnlist[1],expnlist[0],pltvar,sdate,edate,ffmt)
if (fsave):
   print(outname)
   fig.savefig(outname,dpi=quality)
   plt.close()


