#!/usr/bin/env python3
import sys, os, platform
import numpy as np
import netCDF4 as nc
import pandas as pd
import xarray as xa
import matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mpcrs
from datetime import datetime, timedelta
from matplotlib.colors import LinearSegmentedColormap
import setuparea as setarea
from plot_utils import set_size, setupax_2dmap
from utils import find_cnlvs, setup_cmap
import cartopy.crs as ccrs
#
fsave=1 ; ffmt='png' ; ptsize=4
axe_w=6 ; axe_h=3 ; quality=300
tkfreq=3
quality=300; ffmt='png'
cbori='horizontal' #vertical, horizontal
if (cbori=='vertical'):
   cb_frac=0.025
   cb_pad=0.06
elif (cbori=='horizontal'):
   cb_frac=0.04
   cb_pad=0.1

# Projection setting
proj=ccrs.PlateCarree(globe=None)
#
work_year=2016
work_mons=[5,6,7,8,9,10,11,12]
arealist=['Glb','NAmer','Asia','EUR']
hint=6

explist=['noda','anal']
obtype='pm25'

outpath='/data/users/swei/MAPP/'
archdir='/data/users/swei/MAPP/pm25_hofx'
savedir=outpath+'/pm25_2dmap'
suffix='nc4'

if ( not os.path.exists(savedir) ):
    os.makedirs(savedir)

date1 = pd.to_datetime(str(work_year)+'010100',format='%Y%m%d%H')
date2 = pd.to_datetime(str(work_year)+'123118',format='%Y%m%d%H')
delta = timedelta(hours=hint)
dates = pd.date_range(start=date1, end=date2, freq=delta)

nmon=len(work_mons)
if nmon!=1 and nmon<12:
   mons_char=['J','F','M','A','M','J','J','A','S','O','N','D']
   mons_tag=''
   for i,mon_idx in enumerate(work_mons):
       mons_tag+=mons_char[mon_idx-1]
   timetag='%s%s' %(work_year,mons_tag)
elif nmon==12:
   timetag='%s' %(work_year)
else:
   timetag='%s%.2i'%(work_year,work_mons[0])

work_dates=[]
for tmpdate in dates:
    if tmpdate.month in work_mons:
       work_dates.append(tmpdate.strftime('%Y%m%d%H'))
work_dates=pd.to_datetime(work_dates,format='%Y%m%d%H')
ntime=len(work_dates)

def hofx_dict(ncd,varname):
    outdict = {'nlocs':  ncd.dimensions['Location'].size,
               'lats':   ncd.groups['MetaData'].variables['latitude' ][:],
               'lons':   ncd.groups['MetaData'].variables['longitude'][:],
               'obs':    ncd.groups['ObsValue'].variables[varname][:],
               'hofx':   ncd.groups['hofx'].variables[varname][:],
               'preqc':  ncd.groups['PreQC'].variables[varname][:],
               'effqc':  ncd.groups['EffectiveQC'].variables[varname][:],
               }
    return outdict

def process_pm25(in_dict):
    tmpdf = pd.DataFrame.from_dict(in_dict)
    dfout = tmpdf[['effqc','obs','hofx','lats','lons']]
    effqc_filter = (dfout['effqc']==0)
    dfout = dfout.loc[effqc_filter,:]

    return dfout

def calculate_stats(df_in):
    tmp = df_in
    tmp = tmp.reset_index(drop=True)
    tmp['bias'] = tmp['hofx'] - tmp['obs']
    tmp['sq_bias'] = (tmp['hofx'] - tmp['obs'])**2
    
    stats = tmp.groupby(['lats','lons']).agg({'bias':lambda x: x.mean(skipna=True),
                                              'sq_bias': lambda x: x.mean(skipna=True),
                                              'obs': lambda x: x.mean(skipna=True),
                                              'hofx': lambda x: x.mean(skipna=True),
                                             })
    stats['moratio'] = stats['hofx'] / stats['obs']
    stats['corr'] = tmp.groupby(['lats','lons'])[['obs','hofx']].corr().unstack()['obs']['hofx'].values
    stats['corr2'] = stats['corr'].values**2
    stats['rbias'] = stats['bias'] / stats['obs']
    stats['rmse'] = stats['sq_bias'].values**.5
    stats['count'] = tmp.groupby(['lats','lons'])['obs'].count().values

    return stats

for cdate in work_dates:
    cdatestr = cdate.strftime('%Y%m%d%H')
    cymstr = cdate.strftime('%Y%m')

    file1 = ('%s/%s_%s/%s/openaq_pm25_hofx.%s.%s' % (archdir,explist[0],obtype,cymstr,cdatestr,suffix))
    file2 = ('%s/%s_%s/%s/openaq_pm25_hofx.%s.%s' % (archdir,explist[1],obtype,cymstr,cdatestr,suffix))
    if not os.path.exists(file1) or not os.path.exists(file2):
        print('Skipping %s' %cdatestr, flush=1)
        continue

    print('Processing: %s'%file1,flush=1)
    ncd1 = nc.Dataset(file1)
    ncd2 = nc.Dataset(file2)
    dict1 = hofx_dict(ncd1,obtype)
    dict2 = hofx_dict(ncd2,obtype)

    df1 = process_pm25(dict1)
    df2 = process_pm25(dict2)

    if cdate == work_dates[0]:
       df1_all = df1
       df2_all = df2
    else:
       df1_all = pd.concat((df1_all,df1))
       df2_all = pd.concat((df2_all,df2))

stats1 = calculate_stats(df1_all)
stats2 = calculate_stats(df2_all)

stats_name_dict = {'count': 'Counts',
                   'bias': 'Biases',
                   'rmse': 'RMSE',
                   'corr': 'Correlation',
                   'corr2': 'R\u00b2',
                   'moratio': 'Model/Obs.'
                  }

for area in arealist:
    print('Plotting %s' %(area))
    minlon, maxlon, minlat, maxlat, crosszero, cyclic=setarea.setarea(area)
    print(minlat,maxlat,minlon,maxlon,crosszero,cyclic)
    if (area=='Glb'):
       minlon=-180. ; maxlon=180.
    cornll=[minlat,maxlat,minlon,maxlon]
    
    for stat in ['count','bias','rmse','corr','corr2','moratio']:
        if stat == 'bias' or stat == 'corr':
            eqside = 1
            ref_clrmap = 'BlueYellowRed'
        else:
            eqside = 0
            ref_clrmap = 'cmocean_deep'
        cblabel = stats_name_dict[stat]
    
        tmpdf = pd.concat((stats1[stat],stats2[stat]), axis=1)
        tmpdf.columns = explist
        tmpdf.reset_index(inplace=True)
        if area != 'Glb':
            filter = ((tmpdf['lats']<=maxlat)&(tmpdf['lats']>=minlat)&
                      (tmpdf['lons']<=maxlon)&(tmpdf['lons']>=minlon))
            pltdf = tmpdf.loc[filter,:]
        else:
            pltdf = tmpdf
        lats = pltdf['lats'].values
        lons = pltdf['lons'].values
    
        cnlvs = find_cnlvs(pltdf[explist].values,ntcks=21,topq=0.95,eqside=eqside)
        clridx = []
        for idx in np.linspace(2,255,cnlvs.size):
            clridx.append(int(idx))
        clrmap=setup_cmap(ref_clrmap,clridx)
        norm = mpcrs.BoundaryNorm(cnlvs,len(clridx)+1,extend='both')
       
        for exp in explist:
            pltdata = pltdf[exp].values 
            fig,ax,gl=setupax_2dmap(cornll,area,proj,12)
            set_size(axe_w,axe_h)
            sc=ax.scatter(lons,lats,c=pltdata,s=ptsize,cmap=clrmap,norm=norm)
            plt.colorbar(sc,orientation=cbori,label=cblabel,fraction=cb_frac,
                         pad=cb_pad,ticks=cnlvs[::tkfreq],aspect=40)
         
            if fsave:
               outname='%s/%s_%s_%s_%s.%s'%(savedir, area, exp, timetag, stat, ffmt)
               print(outname)
               fig.savefig(outname, dpi=300)
               plt.close()
        
