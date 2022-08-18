__all__ = ['read_rad_ncdiag','read_cnv_ncdiag']

import numpy as np
import xarray as xa

def read_rad_ncdiag(infile,select_wavenumber=None):
    ds=xa.open_dataset(infile)
    npts=int(ds.nobs.size/ds.nchans.size)
    nchs=ds.nchans.size
    rlat=np.reshape(ds.Latitude.values,(npts,nchs))
    rlon=np.reshape(ds.Longitude.values,(npts,nchs))
    rlon=(rlon+180)%360-180
    qcflags=np.reshape(ds.QC_Flag.values,(npts,nchs))
    obs=np.reshape(ds.Observation.values,(npts,nchs))
    sim=np.reshape(ds.Simulated_Tb.values,(npts,nchs))
    clr=np.reshape(ds.Clearsky_Tb.values,(npts,nchs))
    varinv=np.reshape(ds.Inverse_Observation_Error.values,(npts,nchs))
    sim_nbc=np.reshape(ds.Obs_Minus_Forecast_unadjusted.values,(npts,nchs))
    tmpds=xa.Dataset({'rlon':(['obsloc'],rlon[:,0]),
                      'rlat':(['obsloc'],rlat[:,0]),
                      'qcflag':(['obsloc','wavenumber'],qcflags),
                      'tb_obs':(['obsloc','wavenumber'],obs),
                      'tb_sim':(['obsloc','wavenumber'],sim),
                      'tb_clr':(['obsloc','wavenumber'],clr),
                      'varinv':(['obsloc','wavenumber'],varinv),
                      },
                     coords={'obsloc':np.arange(npts),
                             'wavenumber':ds.wavenumber.values})
    if select_wavenumber != None:
       tmpds=tmpds.sel(wavenumber=select_wavenumber)

    return tmpds

def read_cnv_ncdiag(infile):
    ds=xa.open_dataset(infile)
    npts=ds.nobs.size
    rlat=ds.Latitude.data
    rlon=ds.Longitude.data
    rlon=(rlon+180)%360-180
    sta_id=ds.Station_ID.data
    obstype=ds.Observation_Type.data
    sta_elev=ds.Station_Elevation.data
    qcflags=ds.Analysis_Use_Flag.data
    errinv=ds.Errinv_Final.data
    obs=ds.Observation.data
    omb_bc=ds.Obs_Minus_Forecast_adjusted.data
    omb_nbc=ds.Obs_Minus_Forecast_unadjusted.data

    tmpds=xa.Dataset({'rlon':(['obsloc'],rlon),
                      'rlat':(['obsloc'],rlat),
                      'qcflag':(['obsloc'],qcflags),
                      'obs':(['obsloc'],obs),
                      'omb_bc':(['obsloc'],omb_bc),
                      'omb_nbc':(['obsloc'],omb_nbc),
                      'sta_id':(['obsloc'],sta_id),
                      'obstype':(['obsloc'],obstype),
                      },
                     coords={'obsloc':np.arange(npts)})


