#!/usr/bin/env python3

import sys, os
import datetime
import pandas as pd
import pytz as tz
import numpy as np
import netCDF4 as nc
import requests

import lib_python.ioda_conv_engines as iconv
from collections import defaultdict, OrderedDict
from lib_python.orddicts import DefaultOrderedDict

locationKeyList = [
    ("latitude", "float"),
    ("longitude", "float"),
    ("dateTime", "string"),
    ("elevation", "float"),
]

AttrData = {
    'converter': os.path.basename(__file__),
    'nvars': np.int32(1),
}

DimDict = {
}

class improve_measurement(object):
    def __init__(self, in_data, datenow, obsVar):
        self.in_data = in_data
        self.obsVar = obsVar
        self.varDict = defaultdict(lambda: defaultdict(dict))
        self.outdata = defaultdict(lambda: DefaultOrderedDict(OrderedDict))
        self.varAttrs = DefaultOrderedDict(lambda: DefaultOrderedDict(dict))
        self.datenowstr = datenow

        for varname in self.obsVar.keys():
            iodavar = self.obsVar[varname]
            self.varDict[iodavar]['valKey'] = iodavar, iconv.OvalName()
            self.varDict[iodavar]['errKey'] = iodavar, iconv.OerrName()
            self.varDict[iodavar]['qcKey'] = iodavar, iconv.OqcName()
            self.varAttrs[iodavar, iconv.OvalName()]['coordinates'] = 'longitude latitude'
            self.varAttrs[iodavar, iconv.OerrName()]['coordinates'] = 'longitude latitude'
            self.varAttrs[iodavar, iconv.OqcName()]['coordinates'] = 'longitude latitude'
            self.varAttrs[iodavar, iconv.OvalName()]['units'] = 'µg m-3'
            self.varAttrs[iodavar, iconv.OerrName()]['units'] = 'µg m-3'
            self.outdata[self.varDict[iodavar]['valKey']] = np.array(self.in_data[varname],dtype=np.float32)

        self.outdata[('dateTime', 'MetaData')] = np.array(self.in_data['Datetime_UTC'],dtype=object)
        self.outdata[('latitude', 'MetaData')] = np.array(self.in_data['Latitude'],dtype=np.float32)
        self.outdata[('longitude', 'MetaData')] = np.array(self.in_data['Longitude'],dtype=np.float32)
        self.outdata[('deltaTime', 'MetaData')] = np.array(self.in_data['Delta_time'],dtype=np.float32)
        self.outdata[('stationElevation', 'MetaData')] = np.array(self.in_data['Elevation'],dtype=np.float32)
        self.outdata[('stationIdentification', 'MetaData')] = np.array(self.in_data['SiteCode'],dtype=object)

        DimDict['Location'] = len(self.outdata[('dateTime', 'MetaData')])
        AttrData['Location'] = np.int32(DimDict['Location'])
        AttrData['centerdate'] = np.int32(self.datenowstr)
#

def setup_tz_by_states(state):
    pt_states = ['WA','OR','NV','CA']
    mt_states = ['MT','ID','UT','WY','CO','AZ','NM']
    ct_states = ['ND','SD','NE','KS','OK','TX','MN','IA','MO','AR','LA',
                 'WI','IL','KY','TN','MS','AL']
    et_states = ['MI','IN','GA','FL','OH','WV','VA','NC','SC','NY','PA',
                 'NJ','DE','MD','DC','VT','NH','MA','RI','CT','ME']
    at_states = ['AK']
    ht_states = ['HI']
    if state in pt_states:
       tmptz = tz.timezone('US/Pacific')
    elif state in mt_states:
       tmptz = tz.timezone('US/Mountain')
    elif state in ct_states:
       tmptz = tz.timezone('US/Central')
    elif state in et_states:
       tmptz = tz.timezone('US/Eastern')
    elif state in at_states:
       tmptz = tz.timezone('US/Alaska')
    elif state in ht_states:
       tmptz = tz.timezone('US/Hawaii')
    return tmptz

zipcodefile='/data/users/swei/MAPP/uszips.csv'
zipcode_df = pd.read_csv(zipcodefile,dtype={'zip':str,'county_fips':str})
    
outdir = '/data/users/swei/MAPP/improve'
if ( not os.path.exists(outdir) ):
    os.makedirs(outdir)
improvefile = '/data/users/swei/MAPP/improve.xlsx'

obsVar = { 'OCf_Val' : 'concentration_of_organic_carbon_in_air',
           'ECf_Val' : 'concentration_of_elemental_carbon_in_air',
           'SeaSaltf_Val' : 'concentration_of_seasalt_in_air',
           'SO4f_Val' : 'concentration_of_sulfate_in_air',
          }
varDims = { 'concentration_of_organic_carbon_in_air' : ['Location'],
            'concentration_of_elemental_carbon_in_air' : ['Location'],
            'concentration_of_seasalt_in_air' : ['Location'],
            'concentration_of_sulfate_in_air' : ['Location'],
           }

#cycle in 6 hours

sdate=2016010112
edate=2016123112
hint=24

halfday = datetime.timedelta(hours=12)
utc = datetime.timezone.utc

date1 = pd.to_datetime(sdate,format='%Y%m%d%H')
date2 = pd.to_datetime(edate,format='%Y%m%d%H')
delta = datetime.timedelta(hours=hint)
dates = pd.date_range(start=date1, end=date2, freq=delta)

for datenow in dates:
   
    data = pd.read_excel(improvefile,'Data')
    sites = pd.read_excel(improvefile,'Sites',dtype={'County':str})

    data['Datetime'] = pd.to_datetime(data['Date'],format='%m/%d/%Y')+halfday

    filter = ( (data['POC']==1)&(data['Datetime']==datenow) )
    tmpdf = data.loc[filter,:]

    print(datenow.strftime('%Y-%m-%d %H:%M:%S')+' counts: '+str(len(tmpdf)),flush=1)
    if len(tmpdf) == 0:
       continue

    tmpdf = tmpdf[['SiteCode','Datetime','SiteName','Latitude','Longitude','ECf_Val','OCf_Val','SeaSaltf_Val','SO4f_Val']]
    for col in ['ECf_Val','OCf_Val','SeaSaltf_Val','SO4f_Val']:
        cond = (tmpdf[col] != -999.0)
        tmpdf[col] = tmpdf[col].where(cond,np.nan)

    tmpdf['SO4f_Val'] = tmpdf['SO4f_Val'] * 1.375
    tmpdf['OCf_Val'] = tmpdf['OCf_Val'] * 1.8
    
    elev_list = []
    cdate_list = []
    dtime_list = []
    for tmp_stcode in tmpdf['SiteCode']: 
        filter = (sites['Code']==tmp_stcode)
        tmp_st_info_df = sites.loc[filter,:] 
        if tmp_st_info_df['Country'].ravel()[0]=='US':
            if not tmp_st_info_df['County'].isnull().ravel()[0]:
                zipfilter = (zipcode_df['county_fips']==tmp_st_info_df['County'].ravel()[0])
                st_tz = tz.timezone(zipcode_df.loc[zipfilter,:]['timezone'].unique()[0])
            else:
                st_tz = setup_tz_by_states(tmp_st_info_df['State'].ravel()[0])
        elif tmp_st_info_df['Country'].ravel()[0]=='KR':
            st_tz = tz.timezone('Asia/Seoul')
        elif tmp_st_info_df['Country'].ravel()[0]=='CA':
            if tmp_st_info_df['State'].ravel()[0]=='ON':
                st_tz = tz.timezone('Canada/Eastern')
            elif tmp_st_info_df['State'].ravel()[0]=='AB':
                st_tz = tz.timezone('Canada/Mountain')

        center_local_time = datenow.replace(tzinfo=st_tz)
        delta_time_in_sec = center_local_time.utcoffset().total_seconds()
        center_time_in_utc = center_local_time.astimezone(utc)

        elev_list.append(tmp_st_info_df['Elevation'].ravel()[0])
        cdate_list.append(center_time_in_utc.strftime('%Y-%m-%dT%H:%M:%SZ'))
        dtime_list.append(delta_time_in_sec)
     
    tmpdf['Elevation'] = elev_list
    tmpdf['Datetime_UTC'] = cdate_list
    tmpdf['Delta_time'] = dtime_list

    sub_month_dir = datenow.strftime('%Y%m')
    if ( not os.path.exists(outdir+'/'+sub_month_dir) ):
        os.makedirs(outdir+'/'+sub_month_dir)
    outfile = outdir+'/'+sub_month_dir+'/improve.'+datenow.strftime('%Y%m%d00')+'.nc4'

    #sort out obs that have no lat/lon or other unrealistic values
    #write in IODA v3 at datecenter
   
    var = improve_measurement(tmpdf, datenow.strftime('%Y%m%d%H'), obsVar)
 
    # setup the IODA writer
    writer = iconv.IodaWriter(outfile, locationKeyList, DimDict)

    # write everything out
    writer.BuildIoda(var.outdata, varDims, var.varAttrs, AttrData)

    datenow = datenow + datetime.timedelta(hours=hint)
