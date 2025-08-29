#!/usr/bin/env python3
import os
import sys
from openaq import OpenAQ
import time
import requests
import pandas as pd
from datetime import timedelta
from utils import get_dates
import boto3
from httpx import ReadTimeout
from collections import defaultdict

openaq_api_key = '68e3a778022b500e2a88fe70a39bfb33979fcb081b2566788ef8da93af82bea8'
locationstopath = '/work2/noaa/jcsda/shihwei/data/OpenAQ'
if not os.path.exists(locationstopath):
    os.makedirs(locationstopath)

period_sdate = '2018010100'
period_edate = '2022123118'
cycle_interv = 6
window_length = 1
halfwindow = timedelta(hours=window_length/2)

cycle_dates = get_dates(period_sdate, period_edate, cycle_interv)
period_start = pd.to_datetime(period_sdate, format='%Y%m%d%H').tz_localize('UTC') - halfwindow
period_end = pd.to_datetime(period_edate, format='%Y%m%d%H').tz_localize('UTC') + halfwindow

client = OpenAQ(api_key=openaq_api_key)

# Get the latest location list with PM2.5
max_retries = 3
retry_wait = 10  # seconds
locations = []
p = 1
while True:
    print(f'Reading page: {p}', flush=True)
    for attempt in range(max_retries):
        try:
            response = client.locations.list(parameters_id=2, limit=1000, page=p)
            locations += response.results
            break
        except ReadTimeout:
            if attempt < max_retries - 1:
                print(f"ReadTimeout, retrying in {retry_wait} seconds... (attempt {attempt + 1}/{max_retries})", flush=True)
                time.sleep(retry_wait)
            else:
                print("Failed after maximum retries.", flush=True)
                raise

    if response.headers.x_ratelimit_remaining == 0:
        slp_secs = response.headers.x_ratelimit_reset + 1
        print(f"Wait {slp_secs} seconds for limit reset", flush=True)
        time.sleep(slp_secs)

    if response.meta.found == '>1000':
        p += 1
    else:
        break

def keep_available_locations(location, period_start, period_end):
    if location.datetime_first is None or location.datetime_last is None:
        return
    else:
        name = location.name
        id = location.id
        lat = location.coordinates.latitude
        lon = location.coordinates.longitude
        loc_dt_start = pd.to_datetime(location.datetime_first.utc)
        loc_dt_end = pd.to_datetime(location.datetime_last.utc)
        overlapped = (period_end > loc_dt_start) | (loc_dt_end > period_start)
        if overlapped:
            for sensor in location.sensors:
                if 'pm25' in sensor.name:
                    sensor_id = sensor.id
                else:
                    return

    return pd.DataFrame({
        'location': [name],
        'locationID': [id],
        'lat': [lat],
        'lon': [lon],
        'sensorID': [sensor_id],
    })

results = [keep_available_locations(location, period_start, period_end) for location in locations]
locations_df = pd.concat(results).reset_index(drop=True)
print(locations_df.head(), flush=True)

locations_df.to_csv(f'{locationstopath}/tmp.openaq.locations.csv', index=False)
