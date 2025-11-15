#!/usr/bin/env python3

import os
import sys
import glob
import gzip
import pandas as pd
import numpy as np
import netCDF4 as nc
from datetime import datetime, timedelta, timezone
from openaq import OpenAQ

import pyiodaconv.ioda_conv_engines as iconv
import pyiodaconv.ioda_conv_ncio as iconio
from collections import defaultdict, OrderedDict
from pyiodaconv.orddicts import DefaultOrderedDict
from pyiodaconv.def_jedi_utils import iso8601_string, epoch

os.environ["TZ"] = "UTC"

openaq_api_key = '68e3a778022b500e2a88fe70a39bfb33979fcb081b2566788ef8da93af82bea8'
openaq_csv_path = '/work2/noaa/jcsda/shihwei/data/OpenAQ'
availcsv = '20250919.openaq.locations.csv'
outbase = '/work2/noaa/jcsda/shihwei/data/jedi-data/input/obs/openaq_pm25'
period_sdate = '2019010106'
period_edate = '2019013118'
cycle_interv = 6
window_length = 1

obsvars = ['pm25']

varsKeyList = [('valKey', iconv.OvalName(), 'float', 'longitude latitude', 'ug m-3'),
               ('errKey', iconv.OerrName(), 'float', 'longitude latitude', 'ug m-3'),
               ('qcKey', iconv.OqcName(), 'integer', 'longitude latitude', None)]

varDims = {
    'pm25': ['Location'],
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

locationcsv = f'{openaq_csv_path}/{availcsv}'


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


    def create_outdata(self):
        obstime = self.in_data['totalseconds'].values
        obstime = np.where(np.isnan(obstime), long_missing_value, obstime)
        obsval = self.in_data['value'].values
        missing_msk = (np.isnan(obsval) | (obsval < 0.))
        obsval = np.where(missing_msk, float_missing_value, obsval)

        self.outdata[('dateTime', 'MetaData')] = np.array(obstime, dtype=np.int64)
        self.outdata[('latitude', 'MetaData')] = np.array(self.in_data['lat'].values, dtype=np.float32)
        self.outdata[('longitude', 'MetaData')] = np.array(self.in_data['lon'].values, dtype=np.float32)
        self.outdata[('stationIdentification', 'MetaData')] = np.array(self.in_data['location_id'].values.astype(np.str_), dtype=object)
        nloc = len(self.outdata[('dateTime', 'MetaData')])

        for iodavar in obsvars:
            self.outdata[self.varDict[iodavar]['valKey']] = np.array(obsval, dtype=np.float32)
            self.outdata[self.varDict[iodavar]['errKey']] = np.zeros(nloc,dtype=np.float32)
            self.outdata[self.varDict[iodavar]['qcKey']] = np.zeros(nloc,dtype=np.int32)

        DimDict['Location'] = nloc
    
    @staticmethod
    def setStations(stationcsv):
        locations_df = pd.read_csv(stationcsv)
        return locations_df

    @staticmethod
    def load_openaq_data(base_path, dtlist):
        """
        Searches for and loads gzipped CSV files for a specific datetime.
    
        Args:
            base_path (str): The root directory.
            dtlist (datetime.datetime): The specific datetime to search for.
    
        Returns:
            pandas.DataFrame or None: A concatenated DataFrame of all found files,
                                      or None if no files are found.
        """
        file_paths = []
        for dt in dtlist:
            # Extract year, month, and day from the datetime object
            year = dt.year
            month = dt.strftime('%m')  # Format to two digits
            day = dt.strftime('%d')    # Format to two digits
    
            # Construct the base search path
            search_path = os.path.join(base_path, str(year), str(month), '**', f'*-{year}{month}{day}.csv.gz')
            
            # Use glob to find all matching files recursively
            file_paths += glob.glob(search_path, recursive=True)
    
        if not file_paths:
            print(f"No files found for {dt.strftime('%Y-%m-%d')}.", flush=True)
            return None
    
        all_dfs = []
        print(f"Found {len(file_paths)} files. Loading data...", flush=True)
    
        for file_path in file_paths:
            try:
                # Read the gzipped CSV file directly into a DataFrame
                df = pd.read_csv(file_path, compression='gzip')
                all_dfs.append(df)
            except Exception as e:
                print(f"Could not read {file_path}. Error: {e}", flush=True)
    
        # Concatenate all DataFrames into a single one
        if all_dfs:
            return pd.concat(all_dfs, ignore_index=True)
        else:
            return None

def get_month_start_end(year, month):
    # Start of the month
    start = datetime(year, month, 1)
    
    # Find the start of next month
    if month == 12:
        next_month = datetime(year + 1, 1, 1)
    else:
        next_month = datetime(year, month + 1, 1)
    
    # End of this month = 1 second before next month's start
    end = next_month - timedelta(seconds=1)
    
    return start, end

def generate_year_month_list(start_year, start_month, end_year, end_month):
    year_months = []
    year, month = start_year, start_month

    while (year < end_year) or (year == end_year and month <= end_month):
        year_months.append((year, month))
        # Move to next month
        month += 1
        if month > 12:
            month = 1
            year += 1
    return year_months

# 
months=range(1,13)
hours=['00', '06', '12', '18']
client = OpenAQ(api_key=openaq_api_key)
stations = openaq_pm25.setStations(locationcsv)

num_stations = len(stations)
print(num_stations)

for row in stations[['locationID', 'sensorID']].itertuples(index=False):
    climo_holder = df.DataFrame()
    locinfo = client.locations.get(locations_id=row.locationID).results[0]
    loc_first_dt = pd.to_datetime(locinfo.datetime_first.utc)
    loc_last_dt = pd.to_datetime(locinfo.datetime_last.utc)

    year_month_list = generate_year_month_list(
        loc_first_dt.year,
        loc_first_dt.month,
        loc_last_dt.year,
        loc_last_dt.month,
    )
   
    meanval = []
    expected_count = []
    observed_count = []
    avail_yearmonth_list = []
    for year, month in year_month_list:
        dt_from, dt_to = get_month_start_end(year, month)
        measurement = client.measurements.list(
            sensors_id=row.sensorID,
            data="hours",
            rollup="hourofday",
            datetime_from=dt_from,
            datetime_to=dt_to,
            limit=1000,
        )
        if measurement.meta.found != 0:
            avail_yearmonth_list.append((year, month))
            for hourly_mean in measurement.results[0:24:6]:
                meanval.append(hourly_mean.value)
                expected_count.append(hourly_mean.coverage.expected_count)
                observed_count.append(hourly_mean.coverage.observed_count)

        if measurement.headers.x_ratelimit_remaining == 0:
            slp_secs = response.headers.x_ratelimit_reset + 1
            print(f"Wait {slp_secs} seconds for limit reset", flush=True)
            time.sleep(slp_secs)


    expanded = []
    for (y, m) in avail_yearmonth_list:
        for h in hours:
            expanded.append((y, m, h))

    df = pd.DataFrame(expanded, columns=["year", "month", "hour"])
    df['meanval'] = meanval
    df['expected_count'] = expected_count
    df['observed_count'] = observed_count

    station_climatology = df.groupby(["month", "hour"], as_index=False).agg(
        meanval=('meanval', 'mean'),
        expected_count=('expected_count', 'sum'),
        observed_count=('observed_count', 'sum'),
    )
    station_climatology['coverage'] = station_climatology['observed_count'] / station_climatology['expected_count']

sys.exit()

###
period_dt_start = pd.to_datetime(period_sdate, format='%Y%m%d%H')
period_dt_end = pd.to_datetime(period_edate, format='%Y%m%d%H')
halfwindow = timedelta(hours=window_length/2)
cycle_dates = pd.date_range(period_dt_start, period_dt_end, freq=timedelta(hours=cycle_interv))

for cdate, (window_beg, window_end) in zip(cycle_dates, setup_windows(cycle_dates, halfwindow)):
    print(f"Processing: {cdate.strftime('%Y%m%d%H')}", flush=True)
    if window_beg.day != window_end.day:
        date_list = [cdate - timedelta(days=1), cdate]
    else:
        date_list = [cdate]

    df = openaq_pm25.load_openaq_data(openaq_csv_path, date_list)
    df['datetime_no_offset'] = df['datetime'].str[:19]
    df['datetime_utc'] = pd.to_datetime(df['datetime_no_offset'], utc=True)
    parameter_filter = df['parameter'] == 'pm25'
    window_filter = (df['datetime_utc'] >= window_beg) & (df['datetime_utc'] <= window_end)
    tmpdf = df.loc[window_filter & parameter_filter].reset_index(drop=True)
    tmpdf['totalseconds'] = (tmpdf['datetime_utc'] - epoch).dt.total_seconds().astype(np.int64)

    stations = openaq_pm25.setStations(locationcsv)
    stations = stations.rename(columns={'locationID':'location_id'})
    outdf = pd.merge(stations, tmpdf[['location_id', 'parameter', 'units', 'value', 'datetime_utc', 'totalseconds']],
                     on="location_id", how='left')

    outdir = f'{outbase}/{cdate.year}/{cdate.month:02d}'
    if ( not os.path.exists(outdir) ):
        os.makedirs(outdir)

    outfile = outdir+'/openaq_pm25.'+cdate.strftime('%Y%m%d%H')+'.nc4'

    #sort out obs that have no lat/lon or other unrealistic values
    #write in IODA v3 at datecenter
   
    pm25 = openaq_pm25(outdf)
 
    # setup the IODA writer
    writer = iconv.IodaWriter(outfile, locationKeyList, DimDict)

    # write everything out
    writer.BuildIoda(pm25.outdata, varDims, pm25.varAttrs, GlobalAttrs)

