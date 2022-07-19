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
from utils import setup_cmap,find_cnlvs
import cartopy.crs as ccrs

tlsize=10 ; txsize=10
mpl.rc('axes', titlesize=tlsize,labelsize=txsize)
mpl.rc('xtick',labelsize=txsize)
mpl.rc('ytick',labelsize=txsize)
mpl.rc('legend',fontsize='small')
fsave=1 ; ffmt='png' ; ptsize=4
axe_w=6 ; axe_h=3 ; quality=300
tkfreq=2
minussign=u'\u2212'

def plt2d_contourf(fig,ax,pltdata,cnlvs,clrmap,norm,units,cb_dict):
    cn=ax.contourf(pltdata.lon,pltdata.lat,pltdata,levels=cnlvs,
                   cmap=clrmap,norm=norm,extend='both')
    titlestr='Mean=%.2f %s, Max=%.2f %s, Min=%.2f %s' %(pltdata.mean(),units,
                                                        pltdata.max(),units,
                                                        pltdata.min(),units)
    ax.set_title(titlestr,loc='left')
    plt.colorbar(cn,orientation=cb_dict['orient'],fraction=cb_dict['frac'],
                 pad=cb_dict['pad'],aspect=40,label=cb_dict['label'])
    return fig,ax

# Projection setting
proj=ccrs.PlateCarree(globe=None)

cnvtype='t'
units='K'
data_sdate=2020060106
data_edate=2020071018
check_sdate=2020060106
check_edate=2020071018
hint=6
bufrtype='all'
plt_fig=1

degres=2.5

explist=np.array(['hazyda_ctrl','hazyda_aero'])
expnlist=['CTL','AER']

if check_sdate!=data_sdate and check_edate!=data_edate:
   grdtype='mean'
else:
   grdtype='time'

if (cnvtype in ['ps','sst','tcp']):
    is_single_layer=1
elif (cnvtype in ['gps','q','t','uv']):
    is_single_layer=0
    select_pres_lev=850

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
cb_dict={'orient':cbori,'frac':cb_frac,'pad':cb_pad}

inpath=rootarch+'/archive/HazyDA/gridded_diag'
outputpath=rootpath+'/DiagFiles/gridded_cnv'
imgsavpath=outputpath+'/2dmap/omb/'+area+'/'+cnvtype
if ( not os.path.exists(imgsavpath) ):
   os.makedirs(imgsavpath)

grdfile0='%s/%s_%s_%s_%s_omb_%.1fx%.1f.%s.%s_%s.nc' %(inpath,expnlist[0],cnvtype,loop,qcflg,degres,degres,grdtype,data_sdate,data_edate)
grdfile1='%s/%s_%s_%s_%s_omb_%.1fx%.1f.%s.%s_%s.nc' %(inpath,expnlist[1],cnvtype,loop,qcflg,degres,degres,grdtype,data_sdate,data_edate)

ds0=xa.open_dataset(grdfile0)
ds1=xa.open_dataset(grdfile1)

if (not is_single_layer):
   ds0=ds0.sel(lev=select_pres_lev)
   ds1=ds1.sel(lev=select_pres_lev)
   cnvtype_str='%s%s' %(cnvtype,select_pres_lev)

if (bufrtype!='all'):
   select_bufrlist=[]
   for tmpbufrtype in bufrtype:
       if (tmpbufrtype in ds0.obstype.data and
           tmpbufrtype in ds1.obstype.data):
           select_bufrlist.append(tmpbufrtype)
   tmpds0=ds0.sel(obstype=select_bufrlist)
   tmpds1=ds1.sel(obstype=select_bufrlist)
else:
   tmpds0=ds0 ; tmpds1=ds1

if grdtype == 'time':
   select_time_slice=slice(pd.to_datetime(check_sdate,format="%Y%m%d%H"),
                           pd.to_datetime(check_edate,format="%Y%m%d%H"))
   tmpds0=tmpds0.sel(time=select_time_slice)
   tmpds1=tmpds1.sel(time=select_time_slice)
   sum_dims=['obstype','time']
else:
   sum_dims=['obstype']

mean0=(tmpds0['omb_bc_mean']*tmpds0['omb_bc_count']).sum(dim=sum_dims)/(tmpds0['omb_bc_count'].sum(dim=sum_dims))
mean1=(tmpds1['omb_bc_mean']*tmpds1['omb_bc_count']).sum(dim=sum_dims)/(tmpds1['omb_bc_count'].sum(dim=sum_dims))

var0=(tmpds0['omb_bc_var']*tmpds0['omb_bc_count']).sum(dim=sum_dims)/(tmpds0['omb_bc_count'].sum(dim=sum_dims))
var1=(tmpds1['omb_bc_var']*tmpds1['omb_bc_count']).sum(dim=sum_dims)/(tmpds1['omb_bc_count'].sum(dim=sum_dims))

rmse0=np.sqrt(var0+mean0*mean0)
rmse1=np.sqrt(var1+mean1*mean1)
#
#
#
if (plt_fig):
   rmse_cnlvs=find_cnlvs(xa.concat((rmse0,rmse1),dim='exps'),ntcks=21,eqside=0)
   clridx=[]
   for idx in np.linspace(2,254,rmse_cnlvs.size):
       clridx.append(int(idx))
   rmse_clmap=setup_cmap('WhiteYellowOrangeRed',clridx)
   rmse_norm = mpcrs.BoundaryNorm(rmse_cnlvs,len(clridx)+1,extend='both')
   #
   rmse_diff_cnlvs=find_cnlvs(rmse1-rmse0,ntcks=21,eqside=1)
   clridx=[]
   for idx in np.linspace(2,254,rmse_diff_cnlvs.size):
       clridx.append(int(idx))
   rmse_diff_clmap=setup_cmap('BlueYellowRed',clridx)
   rmse_diff_norm = mpcrs.BoundaryNorm(rmse_diff_cnlvs,len(clridx)+1,extend='both')
   #
   mean_cnlvs=find_cnlvs(xa.concat((mean0,mean1),dim='exps'),ntcks=21,eqside=1)
   clridx=[]
   for idx in np.linspace(2,254,mean_cnlvs.size):
       clridx.append(int(idx))
   mean_clmap=setup_cmap('BlueYellowRed',clridx)
   mean_norm = mpcrs.BoundaryNorm(mean_cnlvs,len(clridx)+1,extend='both')
   #
   mean_diff_cnlvs=find_cnlvs(mean1-mean0,ntcks=21,eqside=1)
   clridx=[]
   for idx in np.linspace(2,254,mean_diff_cnlvs.size):
       clridx.append(int(idx))
   mean_diff_clmap=setup_cmap('BlueYellowRed',clridx)
   mean_diff_norm = mpcrs.BoundaryNorm(mean_diff_cnlvs,len(clridx)+1,extend='both')
   #
   #
   #
   if grdtype == 'time':
      sdate_fstr=check_sdate
      edate_fstr=check_edate
   else:
      sdate_fstr=data_sdate
      edate_fstr=data_edate
   #
   fig,ax,gl=setupax_2dmap(cornerll,area,proj,lbsize=txsize)
   set_size(axe_w,axe_h,b=0.15,l=0.05,r=0.95)
   cblabel='%s RMS %s [%s]' %(expnlist[0],tlstr,units)
   cb_dict['label']=cblabel
   plt2d_contourf(fig,ax,rmse0,rmse_cnlvs,rmse_clmap,rmse_norm,units,cb_dict)
   outname='%s/%s_%s_RMS_%s.%s_%s.%s' %(imgsavpath,expnlist[0],cnvtype_str,tlstr,sdate_fstr,edate_fstr,ffmt)
   if (fsave): print(outname,flush=1); fig.savefig(outname,dpi=quality) ; plt.close()
   #
   fig,ax,gl=setupax_2dmap(cornerll,area,proj,lbsize=txsize)
   set_size(axe_w,axe_h,b=0.15,l=0.05,r=0.95)
   cblabel='%s RMS %s [%s]' %(expnlist[1],tlstr,units)
   cb_dict['label']=cblabel
   plt2d_contourf(fig,ax,rmse1,rmse_cnlvs,rmse_clmap,rmse_norm,units,cb_dict)
   outname='%s/%s_%s_RMS_%s.%s_%s.%s' %(imgsavpath,expnlist[1],cnvtype_str,tlstr,sdate_fstr,edate_fstr,ffmt)
   if (fsave): print(outname,flush=1); fig.savefig(outname,dpi=quality) ; plt.close()
   #
   fig,ax,gl=setupax_2dmap(cornerll,area,proj,lbsize=txsize)
   set_size(axe_w,axe_h,b=0.15,l=0.05,r=0.95)
   cblabel='%s%s%s RMS %s [%s]' %(expnlist[1],minussign,expnlist[0],tlstr,units)
   cb_dict['label']=cblabel
   plt2d_contourf(fig,ax,rmse1-rmse0,rmse_diff_cnlvs,rmse_diff_clmap,rmse_diff_norm,units,cb_dict)
   outname='%s/%s-%s_%s_RMS_%s.%s_%s.%s' %(imgsavpath,expnlist[1],expnlist[0],cnvtype_str,tlstr,sdate_fstr,edate_fstr,ffmt)
   if (fsave): print(outname,flush=1); fig.savefig(outname,dpi=quality) ; plt.close()
   #
   fig,ax,gl=setupax_2dmap(cornerll,area,proj,lbsize=txsize)
   set_size(axe_w,axe_h,b=0.15,l=0.05,r=0.95)
   cblabel='%s Mean %s [%s]' %(expnlist[0],tlstr,units)
   cb_dict['label']=cblabel
   plt2d_contourf(fig,ax,mean0,mean_cnlvs,mean_clmap,mean_norm,units,cb_dict)
   outname='%s/%s_%s_Mean_%s.%s_%s.%s' %(imgsavpath,expnlist[0],cnvtype_str,tlstr,sdate_fstr,edate_fstr,ffmt)
   if (fsave): print(outname,flush=1); fig.savefig(outname,dpi=quality) ; plt.close()
   #
   fig,ax,gl=setupax_2dmap(cornerll,area,proj,lbsize=txsize)
   set_size(axe_w,axe_h,b=0.15,l=0.05,r=0.95)
   cblabel='%s Mean %s [%s]' %(expnlist[1],tlstr,units)
   cb_dict['label']=cblabel
   plt2d_contourf(fig,ax,mean1,mean_cnlvs,mean_clmap,mean_norm,units,cb_dict)
   outname='%s/%s_%s_Mean_%s.%s_%s.%s' %(imgsavpath,expnlist[1],cnvtype_str,tlstr,sdate_fstr,edate_fstr,ffmt)
   if (fsave): print(outname,flush=1); fig.savefig(outname,dpi=quality) ; plt.close()
   #
   fig,ax,gl=setupax_2dmap(cornerll,area,proj,lbsize=txsize)
   set_size(axe_w,axe_h,b=0.15,l=0.05,r=0.95)
   cblabel='%s%s%s Mean %s [%s]' %(expnlist[1],minussign,expnlist[0],tlstr,units)
   cb_dict['label']=cblabel
   plt2d_contourf(fig,ax,mean1-mean0,mean_diff_cnlvs,mean_diff_clmap,mean_diff_norm,units,cb_dict)
   outname='%s/%s-%s_%s_Mean_%s.%s_%s.%s' %(imgsavpath,expnlist[1],expnlist[0],cnvtype_str,tlstr,sdate_fstr,edate_fstr,ffmt)
   if (fsave): print(outname,flush=1); fig.savefig(outname,dpi=quality) ; plt.close()

