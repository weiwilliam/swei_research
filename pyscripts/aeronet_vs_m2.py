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
sys.path.append(rootgit+'/pyscripts/functions')
from utils import setup_cmap, ndate, latlon_news
from plot_utils import set_size
import xarray as xa
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.markers as mrk
import matplotlib as mpl
from datetime import datetime, timedelta

# Env and Plot setting
mpl.rc('axes',titlesize=18,labelsize=12)
mpl.rc('xtick',labelsize=12)
mpl.rc('ytick',labelsize=12)
mpl.rc('legend',fontsize='large')
axe_w=8; axe_h=3
quality=300
nofill_dtriangle=mrk.MarkerStyle(marker='v',fillstyle='none')


inputpath='/glade/work/dfgrogan/UFS/WM_DTAER/AER'
outputpath=rootpath+'/Dataset/M2vsAERONET/Timeseries2'
if ( not os.path.exists(outputpath) ):
    os.makedirs(outputpath)

#station_list=['PNNL','Fort_McMurray','']

#station='Boulder'

sdate=2020082200
edate=2020093018
hint=6
m2tag='inst3_2d_gas_Nx'
do_mean=1

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

aeronet_arch='/glade/work/swei/common/AOD/AOD20/TMP'
#aeronet_arch='/glade/work/swei/common/AOD/AOD20/ALL_POINTS'
#infile=aeronet_arch+'/19930101_20211106_'+station+'.lev20'
filelist=os.listdir(aeronet_arch)
filelist.sort()

for infile in filelist:
    df=pd.read_csv(aeronet_arch+'/'+infile,skiprows=6,encoding_errors='ignore')
    df_datetime=pd.to_datetime(df['Date(dd:mm:yyyy)']+' '+df['Time(hh:mm:ss)'],format="%d:%m:%Y %H:%M:%S")
    
    tmpdf=df[["AERONET_Site_Name","Site_Latitude(Degrees)","Site_Longitude(Degrees)","AOD_500nm"]]
    station_name=tmpdf["AERONET_Site_Name"].values[0]
    tmpdf.insert(0,"Datetime",df_datetime)
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


