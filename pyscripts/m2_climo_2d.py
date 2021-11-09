import sys
sys.path.append('/Users/weiwilliam/AlbanyWork/Utility/Python3/functions')
from utils import setup_cmap
from plot_utils import setupax_2dmap,set_size
import setuparea as setarea
import xarray as xa
import numpy as np
from datetime import datetime, timedelta
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mpcrs
import cartopy.crs as ccrs

mpl.rc('axes',titlesize=18,labelsize=16)
mpl.rc('xtick',labelsize=16)
mpl.rc('ytick',labelsize=16)
mpl.rc('legend',fontsize='xx-large')
axe_w=10; axe_h=5

proj=ccrs.PlateCarree()
grav=9.80665e+0

inputpath='/Volumes/WD2TB/ResearchData/MERRA2/climo/'
outputpath='/Users/weiwilliam/AlbanyWork/Dataset/MERRA-2/2dmap'

pltvar='dust'

#mmlst=['06','07','08']
#mmlst=['01','02','03','04','05','06','07','08','09','10','11','12']
mmlst=['06','07']
nmonths=len(mmlst)

pltall=0 # 0: total only 1: sub species included

area='Glb'
minlon, maxlon, minlat, maxlat, crosszero, cyclic=setarea.setarea(area)
print(minlat,maxlat,minlon,maxlon,crosszero,cyclic)
if (area=='Glb'):
   minlon=-180. ; maxlon=180.

#minlat=cornerll[0]; maxlat=cornerll[1]
#minlon=cornerll[2]; maxlon=cornerll[3]
#cornerll=[0.,40.,-80.,20.]
cornerll=[minlat,maxlat,minlon,maxlon]

if (pltvar=='dust'):
   varlst=['DU001','DU002','DU003','DU004','DU005']
elif (pltvar=='carbon'):
   varlst=['OCPHILIC','OCPHOBIC','BCPHILIC','BCPHOBIC'] 
nvars=len(varlst)

for mm in mmlst:
    filepath=inputpath+'/merra2.aerclim.2003-2014.m'+mm+'.nc'
    ds=xa.open_dataset(filepath)
    ds=ds.sel(time=ds.time[0])
   
    delp=ds.DELP
    kgkg_kgm2=delp/grav

    for var in varlst:
        tmp=ds[var]*kgkg_kgm2

        if (var==varlst[0]):
           conc=tmp
           total=tmp
        else:
           conc=xa.concat((conc,tmp),dim='bins')
           total=total+tmp

    conc=xa.concat((conc,total),dim='bins')

    if (mm==mmlst[0]):
       totalconc=conc
    else:
       totalconc=xa.concat((totalconc,conc),dim='time')
    
    del(conc)

# Get the partition of bins
cmass=np.sum(totalconc,axis=2)
# 
clridx=[0,11,20,29,38,47,56,65,74,83,92,101,110,119,128]
cn_cmap=setup_cmap('MPL_YlOrBr',clridx)
#fvmin=0.; fvmax=1e-3
cnlvs=[0      ,  2.5e-5,   5e-5, 7.5e-5,
       1.25e-4,  1.5e-4,1.75e-4,   2e-4,
          3e-4,    4e-4,   5e-4,   6e-4,
          8e-4,    1e-3,   3e-3,   5e-3]
#cnlvs=[0      ,  2.5e-6,   5e-6, 7.5e-6,
#       1.25e-5,  1.5e-5,1.75e-5,   2e-5,
#          3e-5,    4e-5,   5e-5,   6e-5,
#          8e-5,    1e-4,   3e-4,   5e-4]
norm = mpcrs.BoundaryNorm(cnlvs,len(cnlvs))
tkfreq=2
cblb='Column mass density [$\mathrm{kg m^{-2}}$]'

if (not pltall):
   nplotlist=[nvars]
else:
   nplotlist=np.arange(nvars+1)

for m in np.arange(nmonths):
    for n in nplotlist:
        if (n<nvars):
           title='m%s %s column mass density' %(mmlst[m],varlst[n])
           outname='%s/%s_m%s_%s_cmass.png'  %(outputpath,area,mmlst[m],varlst[n])
        else:
           title='m%s total %s column mass density' %(mmlst[m],pltvar) 
           outname='%s/%s_m%s_total%s_cmass.png' %(outputpath,area,mmlst[m],pltvar)

        pltdata=cmass[m,n,:,:]
        fig,ax=setupax_2dmap(cornerll,area,proj,lbsize=16.)
        set_size(axe_w,axe_h,b=0.13,l=0.05,r=0.95,t=0.95)
        cn=ax.contourf(pltdata.lon,pltdata.lat,pltdata,levels=cnlvs,cmap=cn_cmap,norm=norm)
        ax.set_title(title)
        plt.colorbar(cn,ax=ax,orientation='horizontal',ticks=cnlvs[::tkfreq],
                  format='%.2e',fraction=0.045,aspect=40,pad=0.08,label=cblb)
        fig.savefig(outname,dpi=300)
        plt.close()
 
