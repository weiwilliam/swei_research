#!/usr/bin/env python3
import os
import sys
import pandas as pd
from datetime import datetime, timedelta, timezone
from utils import get_dates
import boto3
from httpx import ReadTimeout
from collections import defaultdict

savetopath = '/work2/noaa/jcsda/shihwei/data/OpenAQ'
stationlist = '20250919.openaq.locations.csv'
period_sdate = '2022030100'
period_edate = '2022073118'
cycle_interv = 6
window_length = 1
halfwindow = timedelta(hours=window_length/2)

cycle_dates = get_dates(period_sdate, period_edate, cycle_interv)
period_start = pd.to_datetime(period_sdate, format='%Y%m%d%H').tz_localize('UTC') - halfwindow
period_end = pd.to_datetime(period_edate, format='%Y%m%d%H').tz_localize('UTC') + halfwindow

locationcsvfile = f'{savetopath}/{stationlist}'
locations_df = pd.read_csv(locationcsvfile)
print(locations_df.head(), flush=True)

## Setup s3 bucket
s3 = boto3.client('s3')
bucket_name = "openaq-data-archive"

groups = defaultdict(list)
for dt in cycle_dates:
    groups[(dt.year, dt.month)].append(dt)

# Loop over sorted month-year pairs
for (year, month) in sorted(groups):
    print(f"Month: {month}, Year: {year}", flush=True)

    for row in locations_df.itertuples():
        prefix = f"records/csv.gz/locationid={row.locationID}/year={year}/month={month:02d}"
        savepath = f'{savetopath}/{year}/{month:02d}/location_{row.locationID}'
        if not os.path.exists(savepath):
            os.makedirs(savepath)
    
        available_files = s3.list_objects_v2(Bucket=bucket_name, Prefix=prefix)
        for obj in available_files.get('Contents', []):
            key = obj['Key']
            file = os.path.basename(key)
            remote_size = obj['Size']
            local_file = f'{savepath}/{file}'
            if os.path.exists(local_file):
                local_size = os.path.getsize(local_file)
                if local_size == remote_size:
                    print(f"{local_file} exists and size matches ({local_size} bytes). Skipping download.")
                    continue
                else:
                    print(f"{local_file} exists but size differs (local={local_size}, remote={remote_size}). Re-downloading...")

            s3.download_file(bucket_name, key, local_file)
            print(f"[{datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')}] Downloaded:", local_file, flush=True)
        
