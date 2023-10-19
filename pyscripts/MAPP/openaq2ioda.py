#!/usr/bin/env python3

import os
import datetime
import pandas as pd
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
]

AttrData = {
    'converter': os.path.basename(__file__),
    'nvars': np.int32(1),
}

DimDict = {
}

class openaq_pm25(object):
    def __init__(self, in_data, varname, obsVar):
        self.in_data = in_data
        self.obsVar = obsVar
        self.varDict = defaultdict(lambda: defaultdict(dict))
        self.outdata = defaultdict(lambda: DefaultOrderedDict(OrderedDict))
        self.varAttrs = DefaultOrderedDict(lambda: DefaultOrderedDict(dict))

        iodavar = self.obsVar[varname]
        self.varDict[iodavar]['valKey'] = iodavar, iconv.OvalName()
        self.varDict[iodavar]['errKey'] = iodavar, iconv.OerrName()
        self.varDict[iodavar]['qcKey'] = iodavar, iconv.OqcName()
        self.varAttrs[iodavar, iconv.OvalName()]['coordinates'] = 'longitude latitude'
        self.varAttrs[iodavar, iconv.OerrName()]['coordinates'] = 'longitude latitude'
        self.varAttrs[iodavar, iconv.OqcName()]['coordinates'] = 'longitude latitude'
        self.varAttrs[iodavar, iconv.OvalName()]['units'] = 'µg m-3'
        self.varAttrs[iodavar, iconv.OerrName()]['units'] = 'µg m-3'

        self.outdata[('dateTime', 'MetaData')] = np.array(self.in_data['date'],dtype=object)
        nloc = len(self.outdata[('dateTime', 'MetaData')])

        self.outdata[('latitude', 'MetaData')] = np.array(self.in_data['latitude'],dtype=np.float32)
        self.outdata[('longitude', 'MetaData')] = np.array(self.in_data['longitude'],dtype=np.float32)
        self.outdata[self.varDict[iodavar]['valKey']] = np.array(self.in_data['value'],dtype=np.float32)
        self.outdata[('stationIdentification', 'MetaData')] = np.array(self.in_data['locationId'],dtype=object)
        self.outdata[self.varDict[iodavar]['errKey']] = np.zeros(nloc,dtype=np.float32)
        self.outdata[self.varDict[iodavar]['qcKey']] = np.zeros(nloc,dtype=np.int32)

        DimDict['Location'] = nloc
        AttrData['Location'] = np.int32(DimDict['Location'])
        AttrData['date_time_string'] = self.in_data['date_ref']
#

base_url='https://api.openaq.org/v2/measurements'

outdir = '/data/users/swei/MAPP/openaq_pm25'
if ( not os.path.exists(outdir) ):
    os.makedirs(outdir)

obsVar = { 'pm25' : 'pm25' }
varDims = { 'pm25' : ['Location'] }

#cycle in 6 hours

year_start=2016
month_start=5
day_start=31
hour_start=23
minute_start=30

datestart=datetime.datetime(year_start,month_start,day_start,
                            hour_start,minute_start)
dateend=datestart+datetime.timedelta(days=3)

datenow=datestart
while datenow <= dateend:

    date_from=datenow.strftime("%Y-%m-%dT%X")
    datenowp=datenow+datetime.timedelta(hours=1)
    date_to=datenowp.strftime("%Y-%m-%dT%X")

    params = {
        'date_from':date_from,
        'date_to':date_to,
        'limit':'10000',
        'parameter':'pm25',
        'order_by':'datetime',
        'sensorType':'reference grade',
        'has_geo':'true',
        'value_from':"0",
    }
    
    response = requests.get(base_url,params=params)
    global_data = response.json()

    in_dict = defaultdict(list)
    print('length: %s' % len(global_data['results']))
    for row in global_data['results']:
        for key in ['coordinates','date','locationId','location','value']:
            if key == 'coordinates': 
                if row[key] != None:
                    in_dict['latitude'].append(row[key]['latitude'])
                    in_dict['longitude'].append(row[key]['longitude'])
                else:
                    break
            elif key == 'date':
                in_dict[key].append(row[key]['utc'][:19]+'Z')
            elif key == 'locationId':
                in_dict[key].append(str(row[key]))
            else:
                in_dict[key].append(row[key])
    
    #print(measurements[10])
    datecenter=datenow+datetime.timedelta(minutes=30)
    print(datecenter)
    in_dict['date_ref'] = datecenter.strftime('%Y-%m-%dT%H:%M:%SZ')

    outfile = outdir+'/openaq_pm25.'+datecenter.strftime('%Y%m%d%H')+'.nc4'

    #sort out obs that have no lat/lon or other unrealistic values
    #write in IODA v3 at datecenter
   
    var = openaq_pm25(in_dict,'pm25',obsVar)
 
    # setup the IODA writer
    writer = iconv.IodaWriter(outfile, locationKeyList, DimDict)

    # write everything out
    writer.BuildIoda(var.outdata, varDims, var.varAttrs, AttrData)

    datenow = datenow + datetime.timedelta(hours=6)
