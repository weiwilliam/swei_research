'''
species aod, cmass, etc.
 https://opendap.nccs.nasa.gov/dods/GEOS-5/fp/0.25_deg/assim/inst1_2d_hwl_Nx

species mixing ratio
 https://opendap.nccs.nasa.gov/dods/GEOS-5/fp/0.25_deg/assim/inst3_3d_aer_Nv

'''
import sys
sys.path.append('/Users/weiwilliam/AlbanyWork/Utility/Python3/functions')
import xarray as xa
import numpy as np
from datetime import datetime, timedelta
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mpcrs
from plot_utils import setupax_2dmap, set_size
from utils import setup_cmap
import setuparea as setarea
import cartopy.crs as ccrs

mpl.rc('axes',titlesize=12,labelsize=12)
mpl.rc('xtick',labelsize=12)
mpl.rc('ytick',labelsize=12)
mpl.rc('legend',fontsize='xx-large')
axe_w=8; axe_h=4 ; quality=300 ; lbsize=12

proj=ccrs.PlateCarree()
outputpath='/Users/weiwilliam/AlbanyWork/Dataset/GEOS-5'

grav=9.80665e+0

clridx=[0,11,20,29,38,47,56,65,74,83,92,101,110,119,128]
cn_cmap=setup_cmap('MPL_YlOrBr',clridx)
cnlvs=[0      ,  2.5e-5,   5e-5, 7.5e-5,
       1.25e-4,  1.5e-4,1.75e-4,   2e-4,
          3e-4,    4e-4,   5e-4,   6e-4,
          8e-4,    1e-3,   3e-3,   5e-3]

date=2020091212
aername='carbon'
area='NAmer'
minlon, maxlon, minlat, maxlat, crosszero, cyclic=setarea.setarea(area)
print(minlat,maxlat,minlon,maxlon,crosszero,cyclic)
if (area=='Glb'):
   minlon=-180. ; maxlon=180.
else:
   minlon=(minlon+180)%360-180
   maxlon=(maxlon+180)%360-180
cornerll=[minlat,maxlat,minlon,maxlon]

#cornerll=[-10.,40.,-90.,20.]
#minlat=cornerll[0]; maxlat=cornerll[1]
#minlon=cornerll[2]; maxlon=cornerll[3]

if ( 'du' in aername ):
   aodvar='duexttau'
   varlst=['du001','du002','du003','du004','du005']
if ( 'ss' in aername ):
   aodvar='ssexttau'
   varlst=['ss001','ss002','ss003','ss004','ss005']
if ( 'oc' in aername ):
   aodvar='ocexttau'
   varlst=['ocphobic','ocphilic']
if ( 'bc' in aername ):
   aodvar='bcexttau'
   varlst=['bcphobic','bcphilic']
if ( 'carbon' in aername ):
   varlst=['ocphobic','ocphilic','bcphobic','bcphilic']
if ( 'so4' in aername ):
   aodvar='suexttau'
   varlst=['so4']
nvars=len(varlst)

yy=str(date)[:4]; mm=str(date)[4:6]
dd=str(date)[6:8]; hh=str(date)[8:10]
nc_date = np.datetime64('%s-%s-%sT%s:00:00'%(yy,mm,dd,hh))

ds=xa.open_dataset('https://opendap.nccs.nasa.gov/dods/GEOS-5/fp/0.25_deg/assim/inst3_3d_aer_Nv')
print('Succeed accessing GEOS-5 dataset')

ds=ds.sel(time=nc_date)

if (area != 'Glb'):
   ds=ds.sel(lat=slice(minlat,maxlat),lon=slice(minlon,maxlon))

print('Cropped area and time')

delp=ds.delp
kgkg_kgm2=delp/grav

for var in varlst:
    tmp=ds[var]*kgkg_kgm2
    if (var==varlst[0]):
       conc=tmp
       total=tmp
    else:
       conc=xa.concat((conc,tmp),dim='bins')
       total=total+tmp

totalconc=xa.concat((conc,total),dim='bins')
#pres=prsl

totalconc=totalconc.assign_coords(bins=np.arange(nvars+1))
#pres=pres.assign_coords(lev=np.arange(1.,ds.lev.size+1))
cmass=np.sum(totalconc,axis=1)
print('Finished concentration conversion and integrated for column mass')

norm = mpcrs.BoundaryNorm(cnlvs,len(cnlvs))
tkfreq=2
cblb='Column mass density [$\mathrm{kg m^{-2}}$]'

for n in np.arange(nvars+1):
    if (n<nvars):
       pltdata=cmass[n,:,:]
       title='%s %s column mass density' %(date,varlst[n])
       outname='%s/2dmap/%s_%s_cmass.png'  %(outputpath,date,varlst[n])
    else:
       pltdata=cmass[n,:,:]
       title='%s total %s column mass density' %(date,aername)
       outname='%s/2dmap/%s_total_%s_cmass.png' %(outputpath,date,aername)

    fig,ax=setupax_2dmap(cornerll,area,proj,lbsize)
    set_size(axe_w,axe_h)
    cn=ax.contourf(pltdata.lon,pltdata.lat,pltdata,levels=cnlvs,cmap=cn_cmap,norm=norm)
    ax.set_title(title)
    plt.colorbar(cn,ax=ax,orientation='horizontal',ticks=cnlvs[::tkfreq],
                 format='%.2e',fraction=0.045,aspect=40,pad=0.08,label=cblb)
    fig.savefig(outname,dpi=quality)
    plt.close()
