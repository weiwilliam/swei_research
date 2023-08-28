#!/usr/bin/env python3
import os,sys
import netCDF4 as nc
import numpy as np

no2_path = '/data/users/swei/r2d2_experiments_s4/996a9e/fb/2021-08-01T03:00:00Z'
co_path = '/data/users/swei/Git/jedi-run/tmpdata/v2'

no2_file = '996a9e.fb.tropomi_s5p_no2_total.2021-08-01T03:00:00Z.PT6H.nc4'
co_file = 'fb.tropomi_s5p_co_total.20210801T06_0000.nc4'

no2_ds = nc.Dataset(no2_path+'/'+no2_file)
co_ds = nc.Dataset(co_path+'/'+co_file)

def get_data(ds,var):
    lat = ds.groups['MetaData'].variables['latitude'][:].ravel()
    lon = ds.groups['MetaData'].variables['longitude'][:].ravel()
    qa = ds.groups['MetaData'].variables['quality_assurance_value'][:].ravel()
    hfx = ds.groups['hofx'].variables[var][:].ravel()
    obs = ds.groups['ObsValue'].variables[var][:].ravel()
    
    nloc = ds.dimensions['Location'].size
    nlev = ds.averaging_kernel_levels
    pres = np.zeros((nloc,nlev),dtype='float')
    ak = np.zeros_like(pres)
    
    for i in np.arange(nlev):
        pres_varname = 'pressure_level_%s' %(i+1)
        pres[:,i] = ds.groups['RtrvlAncData'].variables[pres_varname][:].ravel()
        ak_varname = 'averaging_kernel_level_%s' %(i+1)
        ak[:,i] = ds.groups['RtrvlAncData'].variables[ak_varname][:].ravel()

    dict = { 'lat': lat,
             'lon': lon, 
             'hofx': hfx,
             'obs': obs,
             'pres': pres,
             'ak': ak,
             'qa': qa,
           }    
    return dict

def print_out(in_dict, idx):
    rlat = in_dict['lat'][idx]
    rlon = in_dict['lon'][idx]
    print('At (lat,lon) = (%.2f,%.2f)' %(rlat, rlon))
    print('hofx value = %f, qa_value = %f ' %(in_dict['hofx'][idx], in_dict['qa'][idx]))
    print('ak value = ',in_dict['ak'][idx,:])
    
def find_maxidx(in_dict): 
    maxidx = in_dict['hofx'].argmax()
    return maxidx

co_dict = get_data(co_ds,'carbonmonoxideTotal')
no2_dict = get_data(no2_ds,'nitrogendioxideTotal')

maxidx = co_dict['hofx'].argmax()
co_max_lat = co_dict['lat'][maxidx]
co_max_lon = co_dict['lon'][maxidx]

no2_idx = (abs(no2_dict['lat']-co_max_lat) + abs(no2_dict['lon']-co_max_lon)).argmin()
no2_lat = no2_dict['lat'][no2_idx]
no2_lon = no2_dict['lon'][no2_idx]

#print('co max = %f at (lat,lon)=(%.2f,%.2f)' %(co_dict['hofx'][maxidx], co_max_lat, co_max_lon))
#print(co_dict['pres'][maxidx,:])
#print(co_dict['ak'][maxidx,:])
#print('no2 nearest (lat,lon)=(%.2f,%.2f)' %(no2_lat,no2_lon))
#print(no2_dict['pres'][no2_idx,:])
#print(no2_dict['ak'][no2_idx,:])
