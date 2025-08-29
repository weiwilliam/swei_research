#!/usr/bin/env python3

import os
import glob
import gzip
import pandas as pd
import numpy as np
import netCDF4 as nc
from datetime import datetime, timedelta, timezone

import pyiodaconv.ioda_conv_engines as iconv
import pyiodaconv.ioda_conv_ncio as iconio
from collections import defaultdict, OrderedDict
from pyiodaconv.orddicts import DefaultOrderedDict
from pyiodaconv.def_jedi_utils import iso8601_string, epoch

os.environ["TZ"] = "UTC"

openaq_csv_path = '/work2/noaa/jcsda/shihwei/data/OpenAQ'
outdir = '/work2/noaa/jcsda/shihwei/data/jedi-data/input/obs/openaq_pm25'
period_sdate = '2018010100'
period_edate = '2018013118'
cycle_interv = 6
window_length = 1

obsvars = ['particulatematter2p5Insitu']

varsKeyList = [('valKey', iconv.OvalName(), 'float', 'longitude latitude', 'ug m-3'),
               ('errKey', iconv.OerrName(), 'float', 'longitude latitude', 'ug m-3'),
               ('qcKey', iconv.OqcName(), 'integer', 'longitude latitude', None)]

varDims = {
    'particulatematter2p5Insitu': ['Location'],
}

locationKeyList = [
    ("latitude", "float", "degrees_north"),
    ("longitude", "float", "degrees_east"),
    ("dateTime", "long", iso8601_string),
    ("stationIdentification", "string", ""),
]

# A dictionary of variable dimensions.
DimDict = {}

GlobalAttrs = {
    'converter': os.path.basename(__file__),
    'ioda_version': 2,
    'description': 'OpenAQ data (converted from csv.gz to IODA',
    'source': 'https://openaq-data-archive.s3.amazonaws.com/',
    'sourceFiles': '',
}

metaDataName = iconv.MetaDataName()
obsValName = iconv.OvalName()
obsErrName = iconv.OerrName()
qcName = iconv.OqcName()

float_missing_value = iconv.get_default_fill_val(np.float32)
double_missing_value = iconv.get_default_fill_val(np.float64)
int_missing_value = iconv.get_default_fill_val(np.int32)
long_missing_value = iconv.get_default_fill_val(np.int64)
string_missing_value = iconv.get_default_fill_val(np.str_)

missing_vals = {'string': string_missing_value,
                'integer': int_missing_value,
                'long': long_missing_value,
                'float': float_missing_value,
                'double': double_missing_value}

dtypes = {'string': object,
          'integer': np.int32,
          'long': np.int64,
          'float': np.float32,
          'double': np.float64}


def setup_windows(dates, halfwindow):
    beg = (dates - halfwindow).tz_localize(timezone.utc)
    end = (dates + halfwindow).tz_localize(timezone.utc)
    return zip(beg, end)

class openaq_pm25(object):
    def __init__(self, indataframe):
        self.in_data = indataframe
        self.varDict = defaultdict(lambda: defaultdict(dict))
        self.outdata = defaultdict(lambda: DefaultOrderedDict(OrderedDict))
        self.varAttrs = DefaultOrderedDict(lambda: DefaultOrderedDict(dict))
        self.setDicts()
        self.create_outdata()

    def setDicts(self):
        meta_keys = [m_item[0] for m_item in locationKeyList]
        # Set units of the MetaData variables and all _FillValues.
        for key in meta_keys:
            dtypestr = locationKeyList[meta_keys.index(key)][1]
            if locationKeyList[meta_keys.index(key)][2]:
                self.varAttrs[(key, metaDataName)]['units'] = locationKeyList[meta_keys.index(key)][2]
            self.varAttrs[(key, metaDataName)]['_FillValue'] = missing_vals[dtypestr]

        var_keys = [v_item[0] for v_item in varsKeyList]

        # set up variable names for IODA
        for iodavar in obsvars:
            for key in var_keys:
                varGroupName = varsKeyList[var_keys.index(key)][1]
                dtypestr = varsKeyList[var_keys.index(key)][2]
                coord = varsKeyList[var_keys.index(key)][3]
                self.varDict[iodavar][key] = iodavar, varGroupName
                self.varAttrs[iodavar, varGroupName]['coordinates'] = coord
                self.varAttrs[iodavar, varGroupName]['_FillValue'] = missing_vals[dtypestr]
                if varsKeyList[var_keys.index(key)][4]:
                    self.varAttrs[iodavar, varGroupName]['units'] = varsKeyList[var_keys.index(key)][4]
    
    def load_openaq_data(base_path, dt):
        """
        Searches for and loads gzipped CSV files for a specific datetime.
    
        Args:
            base_path (str): The root directory.
            dt (datetime.datetime): The specific datetime to search for.
    
        Returns:
            pandas.DataFrame or None: A concatenated DataFrame of all found files,
                                      or None if no files are found.
        """
        # Extract year, month, and day from the datetime object
        year = dt.year
        month = dt.strftime('%m')  # Format to two digits
        day = dt.strftime('%d')    # Format to two digits
    
        # Construct the base search path
        search_path = os.path.join(base_path, str(year), str(month), '**', f'*-{year}{month}{day}.csv.gz')
        
        # Use glob to find all matching files recursively
        file_paths = glob.glob(search_path, recursive=True)
    
        if not file_paths:
            print(f"No files found for {dt.strftime('%Y-%m-%d')}.")
            return None
    
        all_dfs = []
        print(f"Found {len(file_paths)} files. Loading data...")
    
        for file_path in file_paths:
            try:
                # Read the gzipped CSV file directly into a DataFrame
                df = pd.read_csv(file_path, compression='gzip')
                all_dfs.append(df)
            except Exception as e:
                print(f"Could not read {file_path}. Error: {e}")
    
        # Concatenate all DataFrames into a single one
        if all_dfs:
            return pd.concat(all_dfs, ignore_index=True)
        else:
            return None

    def create_outdata(self):
        self.outdata[('dateTime', 'MetaData')] = np.array(self.in_data['totalseconds'].values, dtype=np.int64)
        self.outdata[('latitude', 'MetaData')] = np.array(self.in_data['lat'].values, dtype=np.float32)
        self.outdata[('longitude', 'MetaData')] = np.array(self.in_data['lon'].values, dtype=np.float32)
        self.outdata[('stationIdentification', 'MetaData')] = np.array(self.in_data['location_id'].values.astype(np.str_), dtype=object)
        nloc = len(self.outdata[('dateTime', 'MetaData')])

        for iodavar in obsvars:
            self.outdata[self.varDict[iodavar]['valKey']] = np.array(self.in_data['value'].values, dtype=np.float32)
            self.outdata[self.varDict[iodavar]['errKey']] = np.zeros(nloc,dtype=np.float32)
            self.outdata[self.varDict[iodavar]['qcKey']] = np.zeros(nloc,dtype=np.int32)

        DimDict['Location'] = nloc

###
if ( not os.path.exists(outdir) ):
    os.makedirs(outdir)

period_dt_start = pd.to_datetime(period_sdate, format='%Y%m%d%H')
period_dt_end = pd.to_datetime(period_edate, format='%Y%m%d%H')
halfwindow = timedelta(hours=window_length/2)
cycle_dates = pd.date_range(period_dt_start, period_dt_end, freq=timedelta(hours=cycle_interv))

for cdate, (window_beg, window_end) in zip(cycle_dates, setup_windows(cycle_dates, halfwindow)):
    print(f"Processing: {cdate.strftime('%Y%m%d%H')}")
    df = openaq_pm25.load_openaq_data(openaq_csv_path, cdate)
    df['datetime_utc'] = pd.to_datetime(df['datetime'], utc=True)
    parameter_filter = df['parameter'] == 'pm25'
    window_filter = (df['datetime_utc'] >= window_beg) & (df['datetime_utc'] <= window_end)
    tmpdf = df.loc[window_filter & parameter_filter].reset_index(drop=True)
    tmpdf['totalseconds'] = (tmpdf['datetime_utc'] - epoch).dt.total_seconds().astype(np.int64)

    outfile = outdir+'/openaq_pm25.'+cdate.strftime('%Y%m%d%H')+'.nc4'

    #sort out obs that have no lat/lon or other unrealistic values
    #write in IODA v3 at datecenter
   
    pm25 = openaq_pm25(tmpdf)
 
    # setup the IODA writer
    writer = iconv.IodaWriter(outfile, locationKeyList, DimDict)

    # write everything out
    writer.BuildIoda(pm25.outdata, varDims, pm25.varAttrs, GlobalAttrs)

