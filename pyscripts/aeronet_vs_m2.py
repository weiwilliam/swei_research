import sys, os, platform
machine='Cheyenne'
if (machine=='MBP'):
    rootpath='/Users/weiwilliam'
    rootarch='/Volumes/WD2TB/ResearchData'
elif (machine=='Desktop'):
    rootpath='F:\GoogleDrive_NCU\Albany'
    rootarch='F:\ResearchData'
    rootgit='F:\GitHub\swei_research'
elif (machine=='S4'):
    rootpath='/data/users/swei'
    rootarch='/scratch/users/swei/ncdiag'
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
from utils import setup_cmap, ndate, latlon_news
from plot_utils import setupax_2dmap, set_size
import setuparea as setarea
import xarray as xa
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.markers as mrk
import matplotlib as mpl
import cartopy.crs as ccrs
import cartopy.feature as cft
from datetime import datetime, timedelta

# Env and Plot setting
mpl.rc('axes',titlesize=18,labelsize=12)
mpl.rc('xtick',labelsize=12)
mpl.rc('ytick',labelsize=12)
mpl.rc('legend',fontsize='large')
axe_w=8; axe_h=3; ptsize=20
quality=300
nofill_dtriangle=mrk.MarkerStyle(marker='v',fillstyle='none')

# Constant configuration
proj=ccrs.PlateCarree()


# Path setup
inputpath='/glade/work/dfgrogan/UFS/WM_DTAER/AER'
outputpath=rootpath+'/Dataset/M2vsAERONET/Timeseries'
if ( not os.path.exists(outputpath) ):
    os.makedirs(outputpath)
sitelistpath=rootpath+'/Dataset/M2vsAERONET/SiteLists'
if ( not os.path.exists(sitelistpath) ):
    os.makedirs(sitelistpath)

#
sdate=2020082200
edate=2020093018
hint=6
m2tag='inst3_2d_gas_Nx'
do_mean=1
area='CONUS'
pltsitemap=1

syy=int(str(sdate)[:4]); smm=int(str(sdate)[4:6])
sdd=int(str(sdate)[6:8]); shh=int(str(sdate)[8:10])
eyy=int(str(edate)[:4]); emm=int(str(edate)[4:6])
edd=int(str(edate)[6:8]); ehh=int(str(edate)[8:10])

date1 = datetime(syy,smm,sdd,shh)
date2 = datetime(eyy,emm,edd,ehh)
delta = timedelta(hours=hint)
delta_15m = timedelta(minutes=15)
dates = pd.date_range(start=date1, end=date2, freq=delta)
date1_m15m=date1-delta_15m
date2_p15m=date2+delta_15m
print('Select AERONET data from '+date1_m15m.strftime('%Y-%m-%d')+' to '+date2_p15m.strftime('%Y-%m-%d'))

tnum=0
dlist=[]
cdate=sdate
while (cdate<=edate):
    dlist.append(str(cdate))
    tnum=tnum+1
    cdate=ndate(hint,cdate)

# Setup plotting area
minlon, maxlon, minlat, maxlat, crosszero, cyclic=setarea.setarea(area)
print(minlat,maxlat,minlon,maxlon,crosszero,cyclic)
if (area=='Glb'):
   minlon=-180. ; maxlon=180.
else:
   minlon=(minlon+180)%360-180
   maxlon=(maxlon+180)%360-180
cornerll=[minlat,maxlat,minlon,maxlon]

#aeronet_arch='/glade/work/swei/common/AOD/AOD20/TMP'
aeronet_arch='/glade/work/swei/common/AOD/AOD20/ALL_POINTS'
#infile=aeronet_arch+'/19930101_20211106_'+station+'.lev20'
filelist=os.listdir(aeronet_arch)
filelist.sort()

prefiltered=0
sitelist_file='%s/%s.%s_%s.txt' %(sitelistpath,area,sdate,edate)
if ( os.path.exists(sitelist_file) ):
   prefilter_df=pd.read_csv(sitelist_file)
   prefiltered=1

staidx=0
for infile in filelist:
    df=pd.read_csv(aeronet_arch+'/'+infile,skiprows=6,encoding_errors='ignore')
    df_datetime=pd.to_datetime(df['Date(dd:mm:yyyy)']+' '+df['Time(hh:mm:ss)'],format="%d:%m:%Y %H:%M:%S")
    
    tmpdf=df[["AERONET_Site_Name","Site_Latitude(Degrees)","Site_Longitude(Degrees)","AOD_500nm"]]
    station_name=tmpdf["AERONET_Site_Name"].values[0]

    if ( prefiltered and not prefilter_df["station"].str.contains(station_name).any() ):
       print('%s is not in the prefiltered list' %(station_name))
       continue

    tmpdf.insert(0,"Datetime",df_datetime)
    if (area!='Glb'):
       filter=((tmpdf["Site_Latitude(Degrees)"]<=maxlat)&(tmpdf["Site_Latitude(Degrees)"]>=minlat)&
               (tmpdf["Site_Longitude(Degrees)"]<=maxlon)&(tmpdf["Site_Longitude(Degrees)"]>=minlon)&
               (tmpdf["Datetime"]<=date2_p15m)&(tmpdf["Datetime"]>=date1_m15m)&
               (tmpdf["AOD_500nm"]>=0.))
    else:
       filter=((tmpdf["Datetime"]<=date2_p15m)&(tmpdf["Datetime"]>=date1_m15m)&
               (tmpdf["AOD_500nm"]>=0.))
    tmpdf=tmpdf.loc[filter,:]
    if (tmpdf.shape[0]==0):
        print('%s : No available data' %(station_name))
        continue

    dates_idx=0
    for date in dlist:
        yy=date[:4] ; mm=date[4:6] ; dd=date[6:8] ; hh=date[8:10]
        pdy=date[:8]
        if (yy=='2020' and mm=='09'):
           m2ind='401'
        else:
           m2ind='400'
    
    # Select AERONET data within +/- 15 minutes
        cdate_p15m=dates[dates_idx]+delta_15m
        cdate_m15m=dates[dates_idx]-delta_15m
        filter=((tmpdf["Datetime"]<=cdate_p15m)&(tmpdf["Datetime"]>=cdate_m15m))
        cdate_df=tmpdf.loc[filter,:]

        if (cdate_df.shape[0]==0):
            stalat=tmpdf["Site_Latitude(Degrees)"].values[0]
            stalon=tmpdf["Site_Longitude(Degrees)"].values[0]
        else: 
            stalat=cdate_df["Site_Latitude(Degrees)"].values[0]
            stalon=cdate_df["Site_Longitude(Degrees)"].values[0]

        if (do_mean):
           cdate_mean_dict={"Datetime":[dates[dates_idx]],"AERONET_Site_Name":[station_name],
                            "Site_Latitude(Degrees)":[stalat],"Site_Longitude(Degrees)":[stalon],
                            "AOD_500nm(Mean)":[cdate_df["AOD_500nm"].mean()]}
           cdate_mean_df=pd.DataFrame(data=cdate_mean_dict)
           if (dates_idx==0):
               aeronet_df=cdate_mean_df
           else:
               aeronet_df=pd.concat((aeronet_df,cdate_mean_df))
        else:    
           if (dates_idx==0):
               aeronet_df=cdate_df
           else:
               aeronet_df=pd.concat((aeronet_df,cdate_df))

    
    # MERRA2_401.inst3_3d_aer_Nv.20200916_12Z.nc4
        infile=inputpath+'/MERRA2_'+m2ind+'.'+m2tag+'.'+pdy+'_'+hh+'Z.nc4'
        ds=xa.open_dataset(infile)
        #tmpaod=ds.AODANA.interp(lat=stalat,lon=stalon,method='linear')
        tmpaod=ds.AODANA.interp(lon=stalon,method='cubic')
        tmpaod=tmpaod.interp(lat=stalat,method='cubic')
    
        if (dates_idx==0):
            m2_aod=tmpaod
        else:
            m2_aod=xa.concat((m2_aod,tmpaod),dim='time')
    
        dates_idx+=1

    if (aeronet_df.shape[0]==0):
        print('%s : No available data within 30 minutes window of each cycle' %(station_name))
        continue
    else:
        print('%s : %i available data' %(station_name,aeronet_df.shape[0]))
    
    txlat,txlon=latlon_news(stalat,stalon)
    titlestr='%s (%s, %s)' %(station_name,txlat,txlon)
    
    fig,ax=plt.subplots()
    set_size(axe_w,axe_h)
    if (do_mean):
       ax.plot_date(aeronet_df["Datetime"],aeronet_df["AOD_500nm(Mean)"],fmt=' ',c='tab:red',marker=nofill_dtriangle)
    else:
       ax.plot_date(aeronet_df["Datetime"],aeronet_df["AOD_500nm"],fmt=' ',c='tab:red',marker=nofill_dtriangle)
    ax.plot_date(dates,m2_aod,fmt='-',c='tab:blue',lw=1.)
    if (do_mean):
       ax.legend(["AOD_500nm(Mean)","M2_AODANA"])
       pltname="AOD_500nm(Mean)"
    else:
       ax.legend(["AOD_500nm","M2_AODANA"])
       pltname="AOD_500nm"
    ax.grid()
    ax.set_ylim(0,5)
    ax.set_title(titlestr,loc='left')
    
    outname='%s/%s_%s.%s_%s.png' %(outputpath,station_name,pltname,sdate,edate)
    
    fig.savefig(outname,dpi=quality)
    plt.close()

    tmp_sites_df=pd.DataFrame(data={'station':station_name,'latitude':stalat,'longitude':stalon},
                              index=[staidx])
    if (staidx==0):
       sites_df=tmp_sites_df
    else:
       sites_df=pd.concat([sites_df,tmp_sites_df])

    if ( not os.path.exists(sitelist_file) ):
       sites_df.to_csv(sitelist_file,index=0)

    staidx+=1

if (pltsitemap):
   fig2,ax2,gl=setupax_2dmap(cornerll,area,proj,lbsize=16.)
   set_size(axe_w,axe_h,b=0.13,l=0.05,r=0.95,t=0.95,ax=ax2)

   ax2.scatter(sites_df['longitude'],sites_df['latitude'],s=ptsize,c='r',marker='s',edgecolors='k')

   for i in np.arange(staidx):
       sta=' %s' %(sites_df['station'].values[i])
       lon=sites_df['longitude'].values[i]
       lat=sites_df['latitude'].values[i]
       ax2.annotate( sta, xy=(lon,lat), zorder=9 )

   gl.xlines=False
   gl.ylines=False
   ax2.add_feature(cft.BORDERS,zorder=2)
   ax2.add_feature(cft.STATES,zorder=2)

   sitemapname='%s/%s_sitemap.%s_%s.png' %(outputpath,area,sdate,edate)
   fig2.savefig(sitemapname,dpi=quality)
   plt.close()
