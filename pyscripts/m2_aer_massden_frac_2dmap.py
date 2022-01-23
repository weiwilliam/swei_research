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
        machine='Cheyenne'
    elif (os.path.exists('/cardinal')):
        rootpath='/data/users/swei/Images'
        rootarch='/scratch2/BMC/gsd-fv3-dev/Shih-wei.Wei/ResearchData'
        rootgit='/home/swei/research'
        machine='S4'
sys.path.append(rootgit+'/pyscripts/functions')
from utils import setup_cmap, ndate
from plot_utils import setupax_2dmap,set_size
import setuparea as setarea
import xarray as xa
import numpy as np
import pandas as pd
from datetime import datetime, timedelta
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mpcrs
import cartopy.crs as ccrs

# Plotting setup
mpl.rc('axes',titlesize=18,labelsize=16)
mpl.rc('xtick',labelsize=16)
mpl.rc('ytick',labelsize=16)
mpl.rc('legend',fontsize='xx-large')
axe_w=10; axe_h=5
quality=300

if (machine=='Cheyenne'):
   inputpath='/glade/work/dfgrogan/UFS/WM_DTAER/AER'
elif (machine=='S4'):
   inputpath='/data/users/swei/common/MERRA2'

outputpath=rootpath+'/Dataset/MERRA-2/2dmap/CMass_frac'
if ( not os.path.exists(outputpath) ):
    os.makedirs(outputpath)

sdate=2020062212
edate=2020062212
hint=6
area='Glb'
m2tag='inst3_3d_aer_Nv'
tkfreq=2
#species_lst=['dust','seas']
species_lst=['dust','seas','carbon','sulf']

# 
clridx=np.array((2,8,9,10,12,14,15,16,18))
cn_cmap=setup_cmap('precip3_16lev',clridx)
cnlvs=np.array((0., 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0))
#clridx=[0,11,29,38,56,74,83,92,119,128]
#cn_cmap=setup_cmap('MPL_YlOrBr',clridx)
#cnlvs=[0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
norm = mpcrs.BoundaryNorm(cnlvs,len(cnlvs))

# Constant configuration
proj=ccrs.PlateCarree()
grav=9.80665e+0

# Setup plotting area
minlon, maxlon, minlat, maxlat, crosszero, cyclic=setarea.setarea(area)
print(minlat,maxlat,minlon,maxlon,crosszero,cyclic)
if (area=='Glb'):
   minlon=-180. ; maxlon=180.
cornerll=[minlat,maxlat,minlon,maxlon]

dust_lst=['DU001','DU002','DU003','DU004','DU005']
seas_lst=['SS001','SS002','SS003','SS004','SS005']
carbon_lst=['OCPHILIC','OCPHOBIC','BCPHILIC','BCPHOBIC'] 
sulf_lst=['SO4']

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

dates_count=0
for date in dlist:
    yy=date[:4] ; mm=date[4:6] ; dd=date[6:8] ; hh=date[8:10]
    pdy=date[:8]
    if (yy=='2020' and mm=='09'):
       m2ind='401'
    else:
       m2ind='400'
# MERRA2_401.inst3_3d_aer_Nv.20200916_12Z.nc4
    if (machine=='Cheyenne'):
       infile=inputpath+'/MERRA2_'+m2ind+'.'+m2tag+'.'+pdy+'_'+hh+'Z.nc4'
       multi_time=0
    if (machine=='S4'):
       infile=inputpath+'/'+yy+'/'+mm+'/MERRA2_'+m2ind+'.'+m2tag+'.'+pdy+'.nc4'
       multi_time=1
   
    ds=xa.open_dataset(infile)
    if (multi_time):
       ds=ds.sel(time=dates[dates_count])
    else:
       ds=ds.sel(time=ds.time[0])
    if (area!='Glb'):
       ds=ds.sel(lon=slice(minlon,maxlon),lat=slice(minlat,maxlat))
   
    delp=ds.DELP
    kgkg_kgm2=delp/grav

    spidx=0
    for sp in species_lst:
        spc_conc=0.
        if (sp=='dust'):
           varlst=dust_lst
        if (sp=='seas'):
           varlst=seas_lst
        if (sp=='carbon'):
           varlst=carbon_lst
        if (sp=='sulf'):
           varlst=sulf_lst
        varidx=0
        for var in varlst:
           tmp=ds[var]*kgkg_kgm2
           if (sp=='dust'):
              if (varidx==0):
                 dust_conc=tmp
              else:
                 dust_conc=xa.concat((dust_conc,tmp),dim='bins')
           if (sp=='seas'):
              if (varidx==0):
                 seas_conc=tmp
              else:
                 seas_conc=xa.concat((seas_conc,tmp),dim='bins')
           if (sp=='carbon'):
              if (varidx==0):
                 carbon_conc=tmp
              else:
                 carbon_conc=xa.concat((carbon_conc,tmp),dim='bins')
           if (sp=='sulf'):
              if (varidx==0):
                 sulf_conc=tmp
              else:
                 sulf_conc=xa.concat((sulf_conc,tmp),dim='bins')
           spc_conc=spc_conc+tmp
           varidx+=1

        if (spidx==0):
           conc=spc_conc
        else:
           conc=xa.concat((conc,spc_conc),dim='species')
        spidx+=1
        
# Get the column integral
    cmass=conc.sum(dim='lev')
    cmass_total=cmass.sum(dim='species')
    
    spidx=0
    for sp in species_lst:
        cblb='%s Fraction of column mass density' %(sp)
        title='%s ' %(species_lst[spidx])
        outname='%s/%s_%s_cmass_frac.%s.png'  %(outputpath,area,species_lst[spidx],date)
    
        pltdata=cmass[spidx,:,:]/cmass_total
        fig,ax=setupax_2dmap(cornerll,area,proj,lbsize=16.)
        set_size(axe_w,axe_h,b=0.13,l=0.05,r=0.95,t=0.95)
        cn=ax.contourf(pltdata.lon,pltdata.lat,pltdata,levels=cnlvs,cmap=cn_cmap,norm=norm)
        ax.set_title(title)
        plt.colorbar(cn,ax=ax,orientation='horizontal',ticks=cnlvs[::tkfreq],
                     format='%.2f',fraction=0.045,aspect=40,pad=0.08,label=cblb)
        print(outname)
        fig.savefig(outname,dpi=quality)
        plt.close()

        if (sp=='dust'):
           varlst=dust_lst
           sp_cmass=dust_conc.sum(dim='lev')
        elif (sp=='seas'):
           varlst=seas_lst
           sp_cmass=seas_conc.sum(dim='lev')
        elif (sp=='carbon'):
           varlst=carbon_lst
           sp_cmass=carbon_conc.sum(dim='lev')
        elif (sp=='sulf'):
           varlst=sulf_lst
           sp_cmass=sulf_conc.sum(dim='lev')
        subspidx=0
        for subsp in varlst:
            cblb='%s Fraction of column mass density' %(subsp)
            title='%s ' %(subsp)
            outname='%s/%s_%s_cmass_frac.%s.png'  %(outputpath,area,subsp,date)
            if (sp=='sulf'):
               pltdata=sp_cmass/cmass_total
            else:
               pltdata=sp_cmass[subspidx,:,:]/cmass_total
            fig,ax=setupax_2dmap(cornerll,area,proj,lbsize=16.)
            set_size(axe_w,axe_h,b=0.13,l=0.05,r=0.95,t=0.95)
            cn=ax.contourf(pltdata.lon,pltdata.lat,pltdata,levels=cnlvs,cmap=cn_cmap,norm=norm)
            ax.set_title(title)
            plt.colorbar(cn,ax=ax,orientation='horizontal',ticks=cnlvs[::tkfreq],
                         format='%.2f',fraction=0.045,aspect=40,pad=0.08,label=cblb)
            print(outname)
            fig.savefig(outname,dpi=quality)
            plt.close()
            subspidx+=1

        spidx+=1
 
    dates_count+=1
