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
from utils import get_dates
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

#
work_year=2016
work_mons=[9,10,11] #3,4,5,6,7,8,9,10,11,12]
arealist=['Glb'] #,'NAmer','Asia','EUR']
hint=6

explist=['noda','anal']
obtype='AOD20'
outpath='/data/users/swei/MAPP/'
aeronetdir='/data/users/swei/MAPP/aeronet_550'
camsdir = '/data/users/swei/MAPP/model/cams'
naradir = '/data/users/swei/MAPP/model/nara'

m2type = 'M2I3NXGAS' # M2T1NXAER
m2_url = 'https://goldsmr4.gesdisc.eosdis.nasa.gov/opendap/MERRA2/%s.5.12.4' %(m2type)
m2tag = 'inst3_2d_gas_Nx' #'tavg1_2d_aer_Nx'
m2ind = 400
savedir=outpath+'/aod550_stats'
wave0 = 550
small = 0.1

if ( not os.path.exists(savedir) ):
    os.makedirs(savedir)

sdate = str(work_year)+'010100'
edate = str(work_year)+'123118'
dates = get_dates(sdate,edate,hint)

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

def hofx_dict(ncd,ds,modelname):
    speed_light = 2.99792458E8 * 1.e9
    wavelengths = ncd.groups['MetaData'].variables['sensorCentralWavelength'][:]
    ichan_wave0 = np.where( abs(wavelengths*1e3 - wave0) < small )[0][0]
    lats = ncd.groups['MetaData'].variables['latitude' ][:]
    lons = ncd.groups['MetaData'].variables['longitude'][:]
    stid = ncd.groups['MetaData'].variables['stationIdentification'][:]
    time = ncd.groups['MetaData'].variables['dateTime'][:].data
    ordinal = np.datetime64('1970-01-01T00:00:00')
    datetime = np.zeros_like(time,dtype='datetime64[s]')
    for n in range(time.size):
        datetime[n] = ordinal + np.timedelta64(time[n],'s')

    st_lats = xa.DataArray(lats,coords=[stid],dims="station")
    st_lons = xa.DataArray(lons,coords=[stid],dims="station")
    st_time = xa.DataArray(datetime,coords=[stid],dims="station")

    # test using select or interp
    if modelname=='nara' or modelname=='merra2':
       hofx = ds.interp(lat=st_lats,lon=st_lons,time=st_time)
    elif modelname=='cams':
       hofx = ds.interp(latitude=st_lats,longitude=st_lons,time=st_time)

    outdict = {
               'lats':   lats,
               'lons':   lons,
               'obs':    ncd.groups['ObsValue'].variables['aerosolOpticalDepth'][:,ichan_wave0],
               'hofx':   hofx,
               }
    return outdict

def process_aod(in_dict):
    tmpdf = pd.DataFrame.from_dict(in_dict)
    dfout = tmpdf[['hofx','obs','lats','lons']]
    nanfilter = (dfout['hofx'].isna()|dfout['obs'].isna())
    dfout = dfout.loc[~nanfilter,:]

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

def datestr(indate):
    datestr = indate.strftime('%Y%m%d%H')
    pdystr = indate.strftime('%Y%m%d')
    ymstr = indate.strftime('%Y%m')
    return datestr, pdystr, ymstr

for d in range(ntime):
    cdate = work_dates[d]
    c_date, c_pdy, c_yms = datestr(cdate)
    print('Processing: %s' %(c_date))
    if d>=1:
        pdate = work_dates[d-1]
        p_date, p_pdy, p_yms = datestr(pdate)
    if d<ntime-1:
        ndate = work_dates[d+1]
        n_date, n_pdy, n_yms = datestr(ndate)
  
    aeronetfile = '%s/aeronet_aod.%s.nc' %(aeronetdir,c_date)
    if not os.path.exists(aeronetfile):
        print('Skipping %s' %c_date, flush=1)
        continue

    aeronet_ncd = nc.Dataset(aeronetfile)
    # CAMS dataset
    if  c_date[-2:]=='00' and d>=1:
        camsfilelist = ['%s/cams_aods_%s.nc' %(camsdir,p_pdy),'%s/cams_aods_%s.nc' %(camsdir,c_pdy)]
    elif c_date[-2:]=='18' and d>=1 and d!=ntime-1:
        camsfilelist = ['%s/cams_aods_%s.nc' %(camsdir,c_pdy),'%s/cams_aods_%s.nc' %(camsdir,n_pdy)]
    else:
        camsfilelist = ['%s/cams_aods_%s.nc' %(camsdir,c_pdy)]
    cams_ds = xa.open_mfdataset(camsfilelist)
    cams_ds = cams_ds['aod550']
    
    # NARA dataset
    if d == 0:
        narafilelist = ['%s/NARA-1.0_AOD_%s.nc4' %(naradir,c_date),'%s/NARA-1.0_AOD_%s.nc4' %(naradir,n_date)]
    elif d == ntime-1:
        narafilelist = ['%s/NARA-1.0_AOD_%s.nc4' %(naradir,p_date),'%s/NARA-1.0_AOD_%s.nc4' %(naradir,c_date)]
    else:
        narafilelist = ['%s/NARA-1.0_AOD_%s.nc4' %(naradir,p_date), \
                    '%s/NARA-1.0_AOD_%s.nc4' %(naradir,c_date), \
                    '%s/NARA-1.0_AOD_%s.nc4' %(naradir,n_date)]
    nara_ds = xa.open_mfdataset(narafilelist)
    nara_ds = nara_ds['AOD'].sel(channel=0)

    # MERRA-2 dataset
    if  c_date[-2:]=='00' and d>=1:
        m2filelist = ['%s/%s/%s/MERRA2_%s.%s.%s.nc4' %(m2_url,p_yms[:4],p_yms[-2:],m2ind,m2tag,p_pdy), \
                      '%s/%s/%s/MERRA2_%s.%s.%s.nc4' %(m2_url,c_yms[:4],c_yms[-2:],m2ind,m2tag,c_pdy), \
                     ]
    elif c_date[-2:]=='18' and d>=1 and d!=ntime-1:
        m2filelist = ['%s/%s/%s/MERRA2_%s.%s.%s.nc4' %(m2_url,c_yms[:4],c_yms[-2:],m2ind,m2tag,c_pdy), \
                      '%s/%s/%s/MERRA2_%s.%s.%s.nc4' %(m2_url,n_yms[:4],n_yms[-2:],m2ind,m2tag,n_pdy), \
                     ]
    else:
        m2filelist = ['%s/%s/%s/MERRA2_%s.%s.%s.nc4' %(m2_url,c_yms[:4],c_yms[-2:],m2ind,m2tag,c_pdy)]
    mer2_ds = xa.open_mfdataset(m2filelist)
    mer2_ds = mer2_ds['AODANA']

    # Dictionary
    cams_dict = hofx_dict(aeronet_ncd,cams_ds,'cams')
    nara_dict = hofx_dict(aeronet_ncd,nara_ds,'nara')
    mer2_dict = hofx_dict(aeronet_ncd,mer2_ds,'merra2')

    cams_df = process_aod(cams_dict)
    nara_df = process_aod(nara_dict)
    mer2_df = process_aod(mer2_dict)

    if d == 0:
       camsdf_all = cams_df
       naradf_all = nara_df
       mer2df_all = mer2_df
    else:
       camsdf_all = pd.concat((camsdf_all,cams_df))
       naradf_all = pd.concat((naradf_all,nara_df))
       mer2df_all = pd.concat((mer2df_all,mer2_df))

cams_stats = calculate_stats(camsdf_all)
nara_stats = calculate_stats(naradf_all)
mer2_stats = calculate_stats(mer2df_all)

allstats = pd.concat((nara_stats.mean(),cams_stats.mean(),mer2_stats.mean()),axis=1)
allstats.columns = ['NARA-1.0','CAMS','MERRA-2']
allstats = allstats.transpose()
outfile = '%s/station_ave.%s.txt' %(savedir,timetag)
allstats.to_string(outfile)

stats_name_dict = {'count': 'Counts',
                   'bias': 'Biases',
                   'rmse': 'RMSE',
                   'corr': 'Correlation',
                   'corr2': 'R\u00b2',
                   'moratio': 'Model/Obs.'
                  }

sys.exit()
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
    
#        cnlvs = find_cnlvs(pltdf[explist].values,ntcks=21,topq=0.95,eqside=eqside)
#        clridx = []
#        for idx in np.linspace(2,255,cnlvs.size):
#            clridx.append(int(idx))
#        clrmap=setup_cmap(ref_clrmap,clridx)
#        norm = mpcrs.BoundaryNorm(cnlvs,len(clridx)+1,extend='both')
#       
#        for exp in explist:
#            pltdata = pltdf[exp].values 
#            fig,ax,gl=setupax_2dmap(cornll,area,proj,12)
#            set_size(axe_w,axe_h)
#            sc=ax.scatter(lons,lats,c=pltdata,s=ptsize,cmap=clrmap,norm=norm)
#            plt.colorbar(sc,orientation=cbori,label=cblabel,fraction=cb_frac,
#                         pad=cb_pad,ticks=cnlvs[::tkfreq],aspect=40)
#         
#            if fsave:
#               outname='%s/%s_%s_%s_%s.%s'%(savedir, area, exp, timetag, stat, ffmt)
#               print(outname)
#               fig.savefig(outname, dpi=300)
#               plt.close()
        
