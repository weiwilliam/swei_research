import os,re
import xarray as xa
import numpy as np

def convert_reflectance(inbandrad,f_sol):
    pi=np.pi
    reflectance=(inbandrad*100*pi)/f_sol
    return reflectance

def convert_brightnessT(chnl_data):
    tb=((c2*wnc)/ln(1.+(c1*wnc^3)/rad)-beta)/alpha
    return tb

def ReadInAVHRR(infile):
    ds=xa.open_dataset(infile,mask_and_scale=1)
    return ds

def ReadInIASIdiag(infile,select_wavenumber=None):
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

def CalAll_ZA(za_first,za,za_last):
    za_tmp=xa.concat((za_first,za),dim='across_track_sampled')
    za_tmp=xa.concat((za_tmp,za_last),dim='across_track_sampled')
    return za_tmp

def main(avhrr_in,ncdiag_in,outimgname,check_wavenumber=None):
    import cartopy.crs as ccrs
    from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
    from pyresample import SwathDefinition
    from pyresample.kd_tree import resample_nearest
    import matplotlib.pyplot as plt
    import matplotlib.colors as mpcrs
    from utils import setup_cmap

    pj_eqc=ccrs.PlateCarree()

    avhrr  = ReadInAVHRR(avhrr_in)
    ncdiag = ReadInIASIdiag(ncdiag_in,check_wavenumber)

    # Get AVHRR channels' data for natural color image
    lat = avhrr.lat.data ; lon=avhrr.lon.data
    minlon=lon.min(); maxlon=lon.max()
    minlat=lat.min(); maxlat=lat.max()
    print(minlon,maxlon,minlat,maxlat)

    solar_za = CalAll_ZA(avhrr.solar_zenith_first,avhrr.solar_zenith,avhrr.solar_zenith_last)
    across_trk_sample_idx=np.zeros(solar_za.shape[1],dtype='int')
    across_trk_sample_idx[1:104]=np.arange(4,2045,20)
    across_trk_sample_idx[104]=2047
    solar_za=solar_za.assign_coords({'across_track_sampled':across_trk_sample_idx,'along_track':avhrr.along_track.data})
    solar_za=solar_za.interp(across_track_sampled=avhrr.across_track.data).rename({'across_track_sampled':'across_track'})
    cos_za=np.cos(np.deg2rad(solar_za))

    ch1_ref  = 0.01*convert_reflectance(avhrr.scene_radiances1,avhrr.channel_1_f_sol) / cos_za
    ch2_ref  = 0.01*convert_reflectance(avhrr.scene_radiances2,avhrr.channel_2_f_sol) / cos_za
    ch3a_ref = 0.01*convert_reflectance(avhrr.scene_radiances3a,avhrr.channel_3a_f_sol) / cos_za
    ch1_ref = ch1_ref.where(ch3a_ref > 0)
    ch2_ref = ch2_ref.where(ch3a_ref > 0)
    ch3a_ref = ch3a_ref.where(ch3a_ref > 0)
    
    rgb_arr = xa.concat((ch3a_ref,ch2_ref,ch1_ref),dim='rgb')
    rgb_arr = rgb_arr.transpose('along_track','across_track','rgb')

    swath_def = SwathDefinition(lon, lat)
    cenlon = 0.5*(minlon+maxlon)
    cenlat = 0.5*(minlat+maxlat)
    area_def = swath_def.compute_optimal_bb_area({"proj": "eqc",'lon_0':cenlon,'lat_0':cenlat})
    crs = area_def.to_cartopy_crs()
    avhrr_rgb = resample_nearest(swath_def, rgb_arr.data, area_def, radius_of_influence=20000, fill_value=None)

    # Define IASI ncdiags
    sc_lon  = ncdiag.rlon ; sc_lat = ncdiag.rlat
    pltmask = (sc_lon<=maxlon)&(sc_lon>=minlon)&(sc_lat<=maxlat)&(sc_lat>=minlat)
    pltlon = sc_lon[pltmask] ; pltlat = sc_lat[pltmask]
    plt_qc= ncdiag.qcflag[pltmask]

    # 
    qc_id_lvs=[0,3,7,10,13,53,55,57]
    qclst=['passed','gross','cloud','phyT','aer_good','sfc_emiss','aer_bad','aer_cld']
    clridx=np.arange(2,2+len(qc_id_lvs))
    clrmap=setup_cmap('grads_default',clridx)
    qc_norm = mpcrs.BoundaryNorm(qc_id_lvs,len(clridx))

    fig = plt.figure()
    ax = plt.subplot(projection=crs)
    ax.set_global()
    ax.coastlines(resolution='110m',color='cyan')
    ax.imshow(avhrr_rgb, extent=crs.bounds, transform=crs, origin='upper')
    sc=ax.scatter(pltlon,pltlat,c=plt_qc,s=3,cmap=clrmap,norm=qc_norm,
                  transform=pj_eqc,zorder=4)
    lg=ax.legend(*sc.legend_elements(),loc=4,title='QC Categories',
                 fontsize=7,title_fontsize=7)

    gl=ax.gridlines(draw_labels=True,dms=True,x_inline=False, y_inline=False)
    gl.right_labels=False
    gl.top_labels=False
    gl.xformatter=LongitudeFormatter(degree_symbol=u'\u00B0 ')
    gl.yformatter=LatitudeFormatter(degree_symbol=u'\u00B0 ')
    gl.xlabel_style={'size':8}
    gl.ylabel_style={'size':8}
    fig.savefig(outimgname,dpi=300)
    plt.close()

if __name__ == "__main__":
    from optparse import OptionParser

    usage = 'usage: %prog -i input-avhrr -n input-ncdiag [-o image_dir -d date_time_group -w assimilation_window]'

    parser = OptionParser(usage)
    parser.add_option('-i', '--input-avhrr', dest='avhrr_file',
                      action='store', type=str, default=None,
                      help='Location of input AVHRR file')
    parser.add_option('-n', '--input-ncdiag', dest='ncdiag_file',
                      action='store', type=str, default=None,
                      help='Location of input nc diag file')
    parser.add_option('-c', '--check-wavenum', dest='check_wavenum',
                      action='store', type=str, default=None,
                      help='wavenumber of target channel')
    parser.add_option('-l', '--exp-label', dest='explabel',
                      action='store', type=str, default=None,
                      help='label of experiment')
    parser.add_option('-d', '--cdtg', dest='cdtg',
                      action='store', type=str, default=None,
                      help='10-digit date time group yyyymmddhh')
    parser.add_option('-o', '--image-dir', dest='image_dir',
                      action='store', type=str, default=os.getcwd(),
                      help='directory path in which to place images')
    (options, args) = parser.parse_args()

    # check for file
    if not options.avhrr_file:
        parser.error("please supply a avhrr file to plot with -i option")
    if not options.ncdiag_file:
        parser.error("please supply a ncdiag file to plot with -n option")
    if not os.path.isfile( options.avhrr_file ):
        print('')
        parser.error("can not find avhrr_file: %s" % options.avhrr_file)
    if not os.path.isfile( options.ncdiag_file ):
        print('')
        parser.error("can not find ncdiag_file: %s" % options.ncdiag_file)

    if not os.path.exists(options.image_dir):
        os.makedirs(options.image_dir)

    fn_lst=re.sub("[/,_]", " ", options.avhrr_file).split()
    for n,tag in enumerate(fn_lst):
       if options.cdtg[:4] in tag:
          break
    time_tag=fn_lst[n]; scan_id=fn_lst[n+1]
    print(time_tag,scan_id)

    imgname='%s_%s_%s.%s_%s.png' %('avhrr_iasi',options.explabel,options.check_wavenum,options.cdtg,scan_id)
    outfilename=os.path.join(options.image_dir,imgname)
    print(outfilename)

    main(options.avhrr_file,options.ncdiag_file,outfilename,
         check_wavenumber=options.check_wavenum)
