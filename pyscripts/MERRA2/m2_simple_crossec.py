#!/usr/bin/env python3

import matplotlib.pyplot as plt
import xarray as xr

m2 = xr.open_dataset('https://goldsmr5.gesdisc.eosdis.nasa.gov/opendap/MERRA2/M2I3NVAER.5.12.4/2019/07/MERRA2_400.inst3_3d_aer_Nv.20190724.nc4')

total_du = m2.DU001 + m2.DU002 + m2.DU003 + m2.DU004 + m2.DU005

region_du = total_du.sel(lon=slice(-60,-20),lat=slice(-5,20))

for i, seltime in enumerate(region_du.time):
    plotdata = region_du.sel(time=seltime).mean(dim='lat')
    fig = plt.figure()
    ax = plt.subplot()
    plotdata.plot(ax=ax)
    ax.invert_yaxis()
    fig.savefig(f'totaldu_{i}.png',dpi=300)
    plt.close()
