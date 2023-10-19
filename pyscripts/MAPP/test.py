#!/usr/bin/env python3
import sys, os, platform
import numpy as np
import netCDF4 as nc
import pandas as pd
import matplotlib
#matplotlib.use('agg')
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mpcrs
from datetime import datetime, timedelta
from matplotlib.colors import LinearSegmentedColormap
from plot_utils import set_size
from utils import find_cnlvs

work_year=2016
work_mons=[1,2,3,4,5,6,7,8,9,10,11,12]
hint=6

explist=['noda','anal']
obtype='AOD20'

outpath='/data/users/swei/MAPP'
archdir='/data/users/swei/MAPP/aeronet_hofx'
savedir=outpath+'/aod_density'

plt_xmax=5.0
plt_in_logbin=1
plt_xybins=30

wave0 = 500
small=0.1
threshold_wave0 = 0.01

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

def hofx_dict(ncd):
    speed_light = 2.99792458E8 * 1.e9
    freqs = ncd.groups['VarMetaData'].variables['frequency'][:]
    outdict = {'nlocs':  ncd.dimensions['nlocs'].size,
               'nchans': ncd.dimensions['nchans'].size,
               'wavelengths': speed_light/freqs,
               'lats':   ncd.groups['MetaData'].variables['latitude' ][:],
               'lons':   ncd.groups['MetaData'].variables['longitude'][:],
               'obs':    ncd.groups['ObsValue'].variables['aerosol_optical_depth'][:],
               'hofx':   ncd.groups['hofx'].variables['aerosol_optical_depth'][:],
               'preqc':  ncd.groups['PreQC'].variables['aerosol_optical_depth'][:],
               'effqc':  ncd.groups['EffectiveQC'].variables['aerosol_optical_depth'][:],
               }
    return outdict

def get_aod(in_dict,wave0):
    ichan_wave0 = np.where( abs(in_dict['wavelengths'].data - wave0) < small )[0][0]
    effqc_mask_wave0 = in_dict['effqc'].data[:,ichan_wave0] == 0
    obs_mask_wave0 = in_dict['obs'].data[:,ichan_wave0] >= threshold_wave0
    effqc_mask = ( effqc_mask_wave0 & obs_mask_wave0 )

    aod_wave0_obs = in_dict['obs'][effqc_mask,ichan_wave0]
    aod_wave0_hofx = in_dict['hofx'][effqc_mask,ichan_wave0]

    return aod_wave0_obs, aod_wave0_hofx, effqc_mask

def calculate_stats(obs,hfx,savedir,expname,timetag):
    correlation_matrix = np.corrcoef(obs, hfx)
    correlation_xy = correlation_matrix[0,1]
    r_squared = correlation_xy**2
    bias=np.mean(hfx)-np.mean(obs)
    moratio=np.mean(hfx)/np.mean(obs)
    rbias=bias/np.mean(obs)
    ssize=len(obs)

    stats_dict = { 'Counts': str("%.0f" % ssize),
                   'Absolute Bias': str("%.3f" % bias),
                   'Relative Bias': str("%.3f" % rbias),
                   'R': str("%.3f" % correlation_xy),
                   'R\u00b2': str("%.3f" % r_squared),
                 }
    print(stats_dict)

    R2str='R\u00b2 = %s' % (str("%.3f" % r_squared))
    biasstr='Bias = %s' % (str("%.3f" % bias))
    sizestr='N = %s' % (str("%.0f" % ssize))

    outfile='%s/Stats_%s_%s_N_Bias_R2_M2ORatio_ReBias.txt'%(savedir,expname,timetag)
    np.savetxt(outfile, [ssize,bias,r_squared, moratio, rbias], delimiter='\b')

    return stats_dict #R2str,biasstr,sizestr

def plt_hist2d(hist2d, axis, stats_dict, xmax, lvs, norm, savedir, expname, timetag):
    xlabstr='AERONET 500 nm AOD'
    ylabstr='GEFS 500 nm AOD'

    fig,ax=plt.subplots()
    set_size(5,5,b=0.1,l=0.1,r=0.95,t=0.95)
    cn=ax.contourf(axis[:-1],axis[:-1],hist2d.swapaxes(0,1), levels=lvs, norm=norm, cmap=white_gist_earth, extend='max')
    plt.plot([0.0, xmax],[0.0, xmax], color='gray', linewidth=2, linestyle='--')
    plt.xlim(0.01, xmax)
    plt.ylim(0.01, xmax)

    if plt_in_logbin:
        ax.set_xscale('log')
        ax.set_yscale('log')
    ax.set_aspect('equal')

    plt.grid(alpha=0.5)
    plt.xlabel(xlabstr, fontsize=11)
    plt.ylabel(ylabstr, fontsize=11)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)

    x_pos = 0.012
    y_pos = 0.988
    for key in stats_dict.keys():
        stat_str = '%s= %s' %(key, stats_dict[key])
        y_pos = y_pos - 0.05
        ax.annotate(stat_str, (x_pos, y_pos), ha='left', va='center', fontsize=12, fontweight='bold', xycoords='axes fraction')

    cb = plt.colorbar(cn,orientation='horizontal',fraction=0.03,aspect=30,pad=0.1, extend='max', ticks=lvs[::50])
    cb.ax.minorticks_off()
    cb.ax.ticklabel_format(axis='x', style='sci', scilimits=(0,0), useMathText=True)
    plt.savefig('%s/%s_%s.png'%(savedir, expname, timetag), format='png',dpi=300)
    plt.close(fig)
    return
    
white_gist_earth = LinearSegmentedColormap.from_list('white_gist_earth', [
    (0,     (1,        1,        1       )),
    (1e-20, (0.965882, 0.915975, 0.913378)),
    (0.2,   (0.772885, 0.646409, 0.444171)),
    (0.4,   (0.568932, 0.677541, 0.340330)),
    (0.6,   (0.249216, 0.576471, 0.342046)),
    (0.8,   (0.143740, 0.396564, 0.488306)),
    (1,     (0.013067, 0.000000, 0.348089)),
    ], N=256)

file2 = ('%s/%s_%s/aeronet_aod_hofx.%s.nc' % (archdir,explist[1],obtype,2016081500))
ncd2 = nc.Dataset(file2)
dict2 = hofx_dict(ncd2)

#for cdate in work_dates:
#    cdatestr = cdate.strftime('%Y%m%d%H')
#
#    file1 = ('%s/%s_%s/aeronet_aod_hofx.%s.nc' % (archdir,explist[0],obtype,cdatestr))
#    file2 = ('%s/%s_%s/aeronet_aod_hofx.%s.nc' % (archdir,explist[1],obtype,cdatestr))
#    if not os.path.exists(file1) or not os.path.exists(file2):
#        print('Skipping %s' %cdatestr, flush=1)
#        continue
#
#    print('Processing: %s'%file2,flush=1)
#    ncd1 = nc.Dataset(file1)
#    ncd2 = nc.Dataset(file2)
#    dict1 = hofx_dict(ncd1)
#    dict2 = hofx_dict(ncd2)
#
#    aod_obs1, aod_hofx1, mask1 = get_aod(dict1,wave0)
#    aod_obs2, aod_hofx2, mask2 = get_aod(dict2,wave0)
#
#    if cdate == work_dates[0]:
#       aod_obs1_all = aod_obs1
#       aod_obs2_all = aod_obs2
#       aod_hofx1_all = aod_hofx1
#       aod_hofx2_all = aod_hofx2
#    else:
#       aod_obs1_all = np.append(aod_obs1_all,aod_obs1)
#       aod_obs2_all = np.append(aod_obs2_all,aod_obs2)
#       aod_hofx1_all = np.append(aod_hofx1_all,aod_hofx1)
#       aod_hofx2_all = np.append(aod_hofx2_all,aod_hofx2)
#
#for ipt in range(2):
#    if ipt == 0:
#        obs=aod_obs1_all
#        hfx=aod_hofx1_all
#    if ipt == 1:
#        obs=aod_obs2_all
#        hfx=aod_hofx2_all
#    expname='%s_%s' % (explist[ipt],obtype)
#    if plt_in_logbin:
#        axis=np.linspace(np.log(0.01),np.log(plt_xmax),plt_xybins+1)
#        axis=np.exp(axis)
#    else:
#        axis=np.linspace(0.01,plt_xmax,plt_xybins+1)
#
#    hist2d,x_edge,y_edge=np.histogram2d(obs, hfx, bins=axis)
#    counts=hist2d.sum()
#    print('%s largest counts: %i' %(timetag,hist2d.max()),flush=1)
#    hist2d/=counts
#    if ipt == 0:
#       data = np.expand_dims(hist2d,axis=0)
#    elif ipt == 1:
#       tmpda = np.expand_dims(hist2d,axis=0)
#       data = np.concatenate((data,tmpda),axis=0)
#
#cnlvs=np.linspace(0,data.max(),256)
#clrnorm = mpcrs.BoundaryNorm(cnlvs,len(cnlvs),extend='max')
#
#for ipt in range(2):
#    if ipt == 0:
#        obs=aod_obs1_all
#        hfx=aod_hofx1_all
#    if ipt == 1:
#        obs=aod_obs2_all
#        hfx=aod_hofx2_all
#    expname='%s_%s' % (explist[ipt],obtype)
#    stats_dict = calculate_stats(obs,hfx,savedir,expname,timetag)
#    plt_hist2d(data[ipt,:,:], axis, stats_dict, plt_xmax, cnlvs, clrnorm, savedir, expname, timetag)
#
