#!/usr/bin/env python3
import sys, os, platform
import numpy as np
import netCDF4 as nc
import pandas as pd
import matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mpcrs
from datetime import datetime, timedelta
from matplotlib.colors import LinearSegmentedColormap
from plot_utils import set_size
#
work_year=2016
work_mons=[1,2,3,4,5,6,7,8,9,10,11,12]
hint=72

explist=['noda','anal']
datatype='improve'
obtype='organic_carbon' # elemental_carbon, organic_carbon, seasalt, sulfate
varlb_dict={'elemental_carbon':'EC',
            'organic_carbon':'OC',
            'seasalt':'Sea Salt',
            'sulfate':'Sulfate'}
varlabel=varlb_dict[obtype]

outpath='/data/users/swei/MAPP/improve_scimg'
archdir='/data/users/swei/MAPP/improve_hofx'
savedir=outpath+'/improve_scatter'

plt_xmax=1
plt_xmin=0.0001
plt_in_logbin=0
plt_xybins=30

def plt_hist2d(hist2d, axis, varlb, stats_dict, xmin, xmax, lvs, norm, savedir, expname, timetag):
    xlabstr='IMPROVE %s ug m-3' %(varlb)
    ylabstr='GEFS %s ug m-3' %(varlb)

    fig,ax=plt.subplots()
    set_size(5,5,b=0.1,l=0.1,r=0.95,t=0.95)
    cn=ax.contourf(axis[:-1],axis[:-1],hist2d.swapaxes(0,1), levels=lvs, norm=norm, cmap=white_gist_earth, extend='max')
    plt.plot([0.0, xmax],[0.0, xmax], color='gray', linewidth=2, linestyle='--')
    plt.xlim(xmin, xmax)
    plt.ylim(xmin, xmax)
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

if ( not os.path.exists(savedir) ):
    os.makedirs(savedir)

date1 = pd.to_datetime(str(work_year)+'010112',format='%Y%m%d%H')
date2 = pd.to_datetime(str(work_year)+'123112',format='%Y%m%d%H')
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

white_gist_earth = LinearSegmentedColormap.from_list('white_gist_earth', [
    (0,     (1,        1,        1       )),
    (1e-20, (0.965882, 0.915975, 0.913378)),
    (0.2,   (0.772885, 0.646409, 0.444171)),
    (0.4,   (0.568932, 0.677541, 0.340330)),
    (0.6,   (0.249216, 0.576471, 0.342046)),
    (0.8,   (0.143740, 0.396564, 0.488306)),
    (1,     (0.013067, 0.000000, 0.348089)),
    ], N=256)

def hofx_dict(infile,varname):
    tmplist = []
    f = open(infile, 'r')
    lines = f.readlines()
    found_var = 0
    for line in lines:
        if varname in line:
            found_var = 1
            continue
        elif line == ' \n' and found_var == 1:
            found_var = 0
            break
        if  found_var == 1: 
            obs = line.split()[0]
            mod = line.split()[1]
            tmplist.append([obs,mod])

    df = pd.DataFrame(tmplist,columns=['obs','hofx'],dtype=np.float32)
    df = df.dropna()
    f.close()

    improve_obs = df['obs'].ravel()[:]
    improve_hofx = df['hofx'].ravel()[:]
    
    return improve_obs, improve_hofx

def calculate_stats(obs,hfx,savedir,expname,timetag):
    correlation_matrix = np.corrcoef(obs, hfx)
    r = correlation_matrix[0,1]
    r_squared = r**2
    bias=np.mean(hfx)-np.mean(obs)
    moratio=np.mean(hfx)/np.mean(obs)
    rbias=bias/np.mean(obs)
    ssize=len(obs)

    stats_dict = { 'Counts': str("%.0f" % ssize),
                   'Absolute Bias': str("%.3f" % bias),
                   'Relative Bias': str("%.3f" % rbias),
                   'R': str("%.3f" % r),
                   'R\u00b2': str("%.3f" % r_squared),
                 }
    print(stats_dict)

    R2str='R\u00b2 = %s' % (str("%.3f" % r_squared))
    Rstr='R\u00b2 = %s' % (str("%.3f" % r))
    biasstr='Bias = %s' % (str("%.3f" % bias))
    sizestr='N = %s' % (str("%.0f" % ssize))

    outfile='%s/Stats_%s_%s_N_Bias_R_R2_M2ORatio_ReBias.txt'%(savedir,expname,timetag)
    np.savetxt(outfile, [ssize,bias,r,r_squared, moratio, rbias], delimiter='\b')

    return stats_dict

for cdate in work_dates:
    cdatestr = cdate.strftime('%Y%m%d%H')
    cymstr = cdate.strftime('%Y%m')

    file1 = ('%s/%s_%s/improve_hofx.%s.txt' % (archdir,explist[0],datatype,cdatestr))
    file2 = ('%s/%s_%s/improve_hofx.%s.txt' % (archdir,explist[1],datatype,cdatestr))
    if not os.path.exists(file1) or not os.path.exists(file2):
        print('Skipping %s' %cdatestr, flush=1)
        continue

    print('Processing: %s'%file1,flush=1)
    improve_obs1, improve_hofx1 = hofx_dict(file1,obtype)
    improve_obs2, improve_hofx2 = hofx_dict(file2,obtype)

    if cdate == work_dates[0]:
       improve_obs1_all = improve_obs1
       improve_obs2_all = improve_obs2
       improve_hofx1_all = improve_hofx1
       improve_hofx2_all = improve_hofx2
    else:
       improve_obs1_all = np.append(improve_obs1_all,improve_obs1)
       improve_obs2_all = np.append(improve_obs2_all,improve_obs2)
       improve_hofx1_all = np.append(improve_hofx1_all,improve_hofx1)
       improve_hofx2_all = np.append(improve_hofx2_all,improve_hofx2)


for ipt in range(2):
    if ipt == 0:
        obs=improve_obs1_all
        hfx=improve_hofx1_all
    if ipt == 1:
        obs=improve_obs2_all
        hfx=improve_hofx2_all
    expname='%s_%s' % (explist[ipt],obtype)
    if plt_in_logbin:
        axis=np.linspace(np.log(plt_xmin),np.log(plt_xmax),plt_xybins+1)
        axis=np.exp(axis)
    else:
        axis=np.linspace(plt_xmin,plt_xmax,plt_xybins+1)
 
    hist2d,x_edge,y_edge=np.histogram2d(obs, hfx, bins=axis)
    counts=hist2d.sum()
    print('%s largest counts: %i' %(timetag,hist2d.max()),flush=1)
    hist2d/=counts
    if ipt == 0:
       data = np.expand_dims(hist2d,axis=0)
    elif ipt == 1:
       tmpda = np.expand_dims(hist2d,axis=0)
       data = np.concatenate((data,tmpda),axis=0)

cnlvs=np.linspace(0,data.max(),256)
clrnorm = mpcrs.BoundaryNorm(cnlvs,len(cnlvs),extend='max')

for ipt in range(2):
    if ipt == 0:
        obs=improve_obs1_all
        hfx=improve_hofx1_all
    if ipt == 1:
        obs=improve_obs2_all
        hfx=improve_hofx2_all
    expname='%s_%s' % (explist[ipt],obtype)
    stats_dict = calculate_stats(obs,hfx,savedir,expname,timetag)
    plt_hist2d(data[ipt,:,:], axis, varlabel, stats_dict, plt_xmin, plt_xmax, cnlvs, clrnorm, savedir, expname, timetag)
