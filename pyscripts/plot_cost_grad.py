import sys, os, platform
os_name=platform.system()
if (os_name=='Darwin'):
    rootpath='/Users/weiwilliam'
    rootarch='/Volumes/WD2TB/ResearchData'
elif (os_name=='Windows'):
    rootpath='F:\GoogleDrive_NCU\Albany'
    rootarch='F:\ResearchData'
    rootgit='F:\GitHub\swei_research'
elif (os_name=='Linux'):
    if (os.path.exists('/scratch1')):
        rootpath='/scratch2/BMC/gsd-fv3-dev/Shih-wei.Wei'
        rootarch='/scratch2/BMC/gsd-fv3-dev/Shih-wei.Wei/ResearchData'
        rootgit='/home/Shih-wei.Wei/research'
    elif (os.path.exists('/glade')):
        rootpath='/glade/work/swei/output/images'
        rootarch='/scratch2/BMC/gsd-fv3-dev/Shih-wei.Wei/ResearchData'
        rootgit='/glade/u/home/swei/research'
    elif (os.path.exists('/cardinal')):
        rootpath='/data/users/swei/Images'
        rootarch='/scratch2/BMC/gsd-fv3-dev/Shih-wei.Wei/ResearchData'
        rootgit='/home/swei/research'
sys.path.append(rootgit+'/pyscripts/functions')
from utils import setup_cmap, ndate
from plot_utils import set_size
import setuparea as setarea
import xarray as xa
import numpy as np
import pandas as pd
from datetime import datetime, timedelta
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mpcrs
# Plotting setup
mpl.rc('axes',titlesize=14,labelsize=12)
mpl.rc('xtick',labelsize=12)
mpl.rc('ytick',labelsize=12)
mpl.rc('legend',fontsize='medium')
axe_w=4; axe_h=4; quality=300; ffmt='png'
tkfreq=1
clrlst=['black','tab:red','tab:blue']
lstylst=['-','-','-']
mrkrlst=[' ',' ',' ']

inputpath='/data/users/swei/archive'
outputpath='/data/users/swei/Images/CostGrad'

sdate=2020060106
edate=2020060812
hint=6
explist=['hazyda_10ctrl']
iter1=50
iter2=150
itermax=iter1+iter2+2

savedir=outputpath+'/'+explist[0]
if ( not os.path.exists(savedir) ):
    os.makedirs(savedir)

#
syy=int(str(sdate)[:4]); smm=int(str(sdate)[4:6])
sdd=int(str(sdate)[6:8]); shh=int(str(sdate)[8:10])
eyy=int(str(edate)[:4]); emm=int(str(edate)[4:6])
edd=int(str(edate)[6:8]); ehh=int(str(edate)[8:10])

date1 = datetime(syy,smm,sdd,shh)
date2 = datetime(eyy,emm,edd,ehh)
delta = timedelta(hours=hint)
dates = pd.date_range(start=date1, end=date2, freq=delta)

tnum=0
dlist=[]
cdate=sdate
while (cdate<=edate):
    dlist.append(str(cdate))
    tnum=tnum+1
    cdate=ndate(hint,cdate)

xaxis=np.arange(0,itermax)

for date in dlist:
    print(date)
    gsistatfile=inputpath+'/'+explist[0]+'/'+date+'/gsistat.gdas.'+date
    f=open(gsistatfile,'r')

    costfunc=[]
    gradient=[]
    for line in f.readlines():
        if ('cost,grad,step,b,step?' in line):
           costfunc.append(line.split()[4])
           gradient.append(line.split()[5])

    cost=np.array(costfunc,dtype='float')
    grad=np.array(gradient,dtype='float')

    cost_fname='%s/costfunc.%s.%s'%(savedir,date,ffmt)
    fig,ax=plt.subplots()
    set_size(axe_w,axe_h,b=0.13)
    ax.plot(xaxis,cost,c='tab:blue')
    ax.legend(explist)
    fig.savefig(cost_fname,dpi=quality)
    plt.close()

    grad_fname='%s/gradient.%s.%s'%(savedir,date,ffmt)
    fig,ax=plt.subplots()
    set_size(axe_w,axe_h,b=0.13)
    ax.plot(xaxis,grad,c='tab:blue')
    ax.legend(explist)
    ax.set_yscale('log')
    fig.savefig(grad_fname,dpi=quality)
    plt.close()
    
    f.close()
