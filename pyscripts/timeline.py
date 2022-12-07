# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 17:17:39 2022

@author: ck102
"""
import numpy as np
import datetime as dt
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.collections import PolyCollection
from plot_utils import set_size

tlsize=12 ; txsize=12
mpl.rc('axes', titlesize=tlsize,labelsize=txsize)
mpl.rc('xtick',labelsize=txsize)
mpl.rc('ytick',labelsize=txsize)
mpl.rc('legend',fontsize='small')


data = [(dt.datetime(2022,10,23,13,20), dt.datetime(2022,10,23,13,40), 'Prep_v1.1'),
        (dt.datetime(2022,10,23,13,40), dt.datetime(2022,10,23,13,46), 'GSI_v1.1'),
        (dt.datetime(2022,10,23,13,46), dt.datetime(2022,10,23,19,38), 'WRF-Chem_v1.1'),
        (dt.datetime(2022,10,23,19,20), dt.datetime(2022,10,23,19,52), 'Prep_v1.1'),
        (dt.datetime(2022,10,23,19,52), dt.datetime(2022,10,23,19,57), 'GSI_v1.1'),
        (dt.datetime(2022,10,23,19,57), dt.datetime(2022,10,23,20,30), 'WRF-Chem_v1.1'),
        (dt.datetime(2022,10,24, 1,20), dt.datetime(2022,10,24, 1,37), 'Prep_v1.1'),
        (dt.datetime(2022,10,24, 1,37), dt.datetime(2022,10,24, 1,42), 'GSI_v1.1'),
        (dt.datetime(2022,10,24, 1,42), dt.datetime(2022,10,24, 2,21), 'WRF-Chem_v1.1'),
        (dt.datetime(2022,10,24, 7,20), dt.datetime(2022,10,24, 7,37), 'Prep_v1.1'),
        (dt.datetime(2022,10,24, 7,37), dt.datetime(2022,10,24, 7,42), 'GSI_v1.1'),
        (dt.datetime(2022,10,24, 7,42), dt.datetime(2022,10,24, 8,19), 'WRF-Chem_v1.1'),
        (dt.datetime(2022,10,23,13,20), dt.datetime(2022,10,23,14,10), 'Prep_v1.0'),
        (dt.datetime(2022,10,23,14,10), dt.datetime(2022,10,23,14,15), 'GSI_v1.0'), 
        (dt.datetime(2022,10,23,14,15), dt.datetime(2022,10,23,22,25), 'WRF-Chem_v1.0'),
        (dt.datetime(2022,10,23,19,20), dt.datetime(2022,10,23,19,46), 'Prep_v1.0'),
        (dt.datetime(2022,10,23,19,46), dt.datetime(2022,10,23,19,51), 'GSI_v1.0'), 
        (dt.datetime(2022,10,23,19,51), dt.datetime(2022,10,23,21,15), 'WRF-Chem_v1.0'),
        (dt.datetime(2022,10,24, 1,20), dt.datetime(2022,10,24, 1,46), 'Prep_v1.0'),
        (dt.datetime(2022,10,24, 1,46), dt.datetime(2022,10,24, 1,51), 'GSI_v1.0'), 
        (dt.datetime(2022,10,24, 1,51), dt.datetime(2022,10,24, 3,15), 'WRF-Chem_v1.0'),
        (dt.datetime(2022,10,24, 7,20), dt.datetime(2022,10,24, 7,46), 'Prep_v1.0'),
        (dt.datetime(2022,10,24, 7,46), dt.datetime(2022,10,24, 7,51), 'GSI_v1.0'), 
        (dt.datetime(2022,10,24, 7,51), dt.datetime(2022,10,24, 9,15), 'WRF-Chem_v1.0'),
        (dt.datetime(2022,10,23, 8, 0), dt.datetime(2022,10,23, 8, 2), 'CHEM_FINN'),
        (dt.datetime(2022,10,23,13, 0), dt.datetime(2022,10,23,13,27), 'CHEM_WACCM'),
        (dt.datetime(2022,10,23, 9,30), dt.datetime(2022,10,23, 9,32), 'GDAS_OBS'),
        (dt.datetime(2022,10,23,15,30), dt.datetime(2022,10,23,15,32), 'GDAS_OBS'),
        (dt.datetime(2022,10,23,21,30), dt.datetime(2022,10,23,21,32), 'GDAS_OBS'),
        (dt.datetime(2022,10,24, 3,30), dt.datetime(2022,10,24, 3,32), 'GDAS_OBS'),
        (dt.datetime(2022,10,23, 9, 0), dt.datetime(2022,10,23, 9, 2), 'GFS_OBS'),
        (dt.datetime(2022,10,23,15, 0), dt.datetime(2022,10,23,15, 2), 'GFS_OBS'),
        (dt.datetime(2022,10,23,21, 0), dt.datetime(2022,10,23,21, 2), 'GFS_OBS'),
        (dt.datetime(2022,10,24, 3, 0), dt.datetime(2022,10,24, 3, 2), 'GFS_OBS'),
        (dt.datetime(2022,10,23,10, 0), dt.datetime(2022,10,23,10,13), 'GFS_GRIB2'),
        (dt.datetime(2022,10,23,16, 0), dt.datetime(2022,10,23,16,16), 'GFS_GRIB2'),
        (dt.datetime(2022,10,23,22, 0), dt.datetime(2022,10,23,22,14), 'GFS_GRIB2'),
        (dt.datetime(2022,10,24, 4, 0), dt.datetime(2022,10,24, 4,15), 'GFS_GRIB2'),
        ]

cats = {
        "WRF-Chem_v1.0" : 1, "GSI_v1.0" : 2, "Prep_v1.0" : 3,
        "WRF-Chem_v1.1" : 4, "GSI_v1.1" : 5, "Prep_v1.1" : 6, 
        "GFS_GRIB2":7, "GFS_OBS":8, "GDAS_OBS": 9,
        "CHEM_WACCM" :10, "CHEM_FINN":11, 
        }
colormapping = {
        "CHEM_FINN" : "C3", "CHEM_WACCM": "C3", 
        "GDAS_OBS": "C4", "GFS_OBS": "C5", "GFS_GRIB2": "C5",
        "Prep_v1.1" : "C0", "GSI_v1.1" : "C1", "WRF-Chem_v1.1" : "C2",
        "Prep_v1.0" : "C0", "GSI_v1.0" : "C1", "WRF-Chem_v1.0" : "C2"}

verts = []
colors = []
for d in data:
    v =  [(mdates.date2num(d[0]), cats[d[2]]-.4),
          (mdates.date2num(d[0]), cats[d[2]]+.4),
          (mdates.date2num(d[1]), cats[d[2]]+.4),
          (mdates.date2num(d[1]), cats[d[2]]-.4),
          (mdates.date2num(d[0]), cats[d[2]]-.4)]
    verts.append(v)
    colors.append(colormapping[d[2]])

bars = PolyCollection(verts, facecolors=colors, alpha=0.75, edgecolors=colors)

dfmt = mdates.DateFormatter('%Y %h %n %d %Hz')

fig, ax = plt.subplots()
set_size(8,4,l=0.15)
ax.add_collection(bars)
ax.autoscale()
mjloc = mdates.HourLocator(byhour=[0,6,12,18,])
miloc = mdates.HourLocator(byhour=[3,9,15,21,])
ax.xaxis.set_major_locator(mjloc)
ax.xaxis.set_minor_locator(miloc)
ax.xaxis.set_major_formatter(dfmt)
ax.grid(axis='x')

ax.set_yticks(np.arange(1,12))
ax.set_yticklabels(list(cats.keys()))
fig.savefig('/data/users/swei/FTPdir/Wx-AQ_timeline.png',dpi=300)
plt.close()

