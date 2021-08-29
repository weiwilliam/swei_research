#!/bin/sh --login
##SBATCH --output=../logs/gsirun.out
##SBATCH --job-name=swei_gsistandalone
###SBATCH --qos=debug
##SBATCH --qos=batch
##SBATCH --time=2:00:00
##SBATCH --nodes=8 --ntasks-per-node=20 --cpus-per-task=1
##SBATCH --account=gsd-fv3-test

#export OMP_NUM_THREADS=1
set -x

#
# Set experiment name and analysis date
machine='s4'
exp="aer_observer"
expid=2 # 1: no aer 2: aer
VERBOSE='.false.'
if_observer=Yes 

#
if [ $machine == 'hera' ]; then
   homedir=/scratch1/BMC/gsd-fv3-dev/Shih-wei.Wei/SingleRadTest
   scrpts_home=/home/swei/research/GSI_exps # Modify on Hera
   aerpath=/scratch1/BMC/gsd-fv3-dev/Shih-wei.Wei/common/MERRA2
   obsarch=$homedir/INPUTS/OBSDATA
   DMPDIR=/scratch1/NCEPDEV/global/glopara/dump
   . /apps/lmod/lmod/init/sh
   module purge
   source /scratch2/BMC/gsd-fv3-dev/Shih-wei.Wei/GSI/modulefiles/modulefile.ProdGSI.hera
   module list
elif [ $machine == 's4' ]; then
   homedir=/data/users/swei/Experiments/${exp}_GSI
   scrpts_home=/home/swei/research/GSI_exps
   aerpath=/data/users/swei/common/MERRA2_L64
   obsarch=/data/prod/glopara/dump
   gesarch=/data/users/swei/common/ICs
   DMPDIR=/scratch/users/swei/runtmp
   module purge
   source /data/users/swei/Git/GSI/modulefiles/modulefile.ProdGSI.s4
   module list
else
   echo 'not supported machine, exit'
   exit 1
fi

if [ ${if_observer} = Yes ] ; then
  nummiter=0
  if_read_obs_save='.true.'
  if_read_obs_skip='.false.'
else
  nummiter=2
  if_read_obs_save='.false.'
  if_read_obs_skip='.false.'
fi

case $expid in
1)
  READEXTAER='.false.'
   MERRA2AER='.false.'
  satinfo=${scrpts_home}/dat/controlrun_satinfo.txt
 anavinfo=${scrpts_home}/dat/anavinfo_controlrun ;;
2)
  READEXTAER='.true.'
   MERRA2AER='.true.'
  satinfo=${scrpts_home}/dat/fv3aerorad_satinfo.txt
 anavinfo=${scrpts_home}/dat/anavinfo_fv3aerorad ;;
esac

# Set the JCAP resolution which you want.
# All resolutions use LEVS=64
export ACASE=C384
anl_res=$(echo $ACASE | cut -c2-)
export JCAP_A=$((anl_res*2-2))
export NLON_A=$((anl_res*4))
export NLAT_A=$((anl_res*2+2))

export BCASE=C768
fg_res=$(echo $BCASE| cut -c2-)
export LONB=$((fg_res*4))  #3072
export LATB=$((fg_res*2))  #1536
export JCAP=$((LATB-2))
export LEVS=64

# Set variables used in script
#   CLEAN up $tmpdir when finished (YES=remove, NO=leave alone)
#   ncp is cp replacement, currently keep as /bin/cp
CDATE=2020062212

# Set runtime and save directories
tmpdir=${DMPDIR}/rundir${expid}
outdir=${homedir}/OUTPUT/$exp/$CDATE
logsdir=$homedir/logs
[[ ! -d $outdir ]] && mkdir -p $outdir
[[ ! -d $logsdir ]] && mkdir -p $logsdir

endianness="Big_Endian"

if [ $machine == 'hera' ]; then
   gsidir="/scratch2/BMC/gsd-fv3-dev/Shih-wei.Wei/GSI"
   fixcrtm="/scratch1/NCEPDEV/da/Michael.Lueken/CRTM_REL-2.2.3/crtm_v2.2.3/fix_update"
   COMdir="/scratch1/BMC/gsd-fv3-dev/Shih-wei.Wei/common"
elif [ $machine == 's4' ]; then
   gsidir="/data/users/swei/Git/GSI"
   fixcrtm="/data/users/swei/libs/crtm_coeff/v2.3.0"
   #NDATE="${PROD_UTIL}/bin/ndate"
   COMdir="/data/users/swei/common"
fi
ndate=${NDATE:-/scratch2/NCEPDEV/nwprod/NCEPLIBS/utils/prod_util.v1.1.0/exec/ndate}
fixgsi="${gsidir}/fix"
gsiexec=${gsidir}/exec/global_gsi.x
CATEXEC=${gsidir}/exec/ncdiag_cat.x
cdump=gdas

UNCOMPRESS=gunzip
CLEAN=NO
APRUN=srun
ncp=/bin/cp
ncpl="ln -fs"

# Given the analysis date, compute the date from which the
# first guess comes.  Extract cycle and set prefix and suffix
# for guess and observation data files
PDY=`echo $CDATE | cut -c1-8`
ayy=`echo $CDATE | cut -c1-4`
amm=`echo $CDATE | cut -c5-6`
add=`echo $CDATE | cut -c7-8`
cyc=`echo $CDATE | cut -c9-10`

gdate=`$ndate -6 $CDATE`
 gPDY=`echo $gdate | cut -c1-8`
 gcyc=`echo $gdate | cut -c9-10`

aerp3date=`$ndate  3 $CDATE`
     p3yy=`echo $aerp3date | cut -c1-4`
     p3mm=`echo $aerp3date | cut -c5-6`
     p3dd=`echo $aerp3date | cut -c7-8`
     p3hh=`echo $aerp3date | cut -c9-10`
aerm3date=`$ndate -3 $CDATE`
     m3yy=`echo $aerm3date | cut -c1-4`
     m3mm=`echo $aerm3date | cut -c5-6`
     m3dd=`echo $aerm3date | cut -c7-8`
     m3hh=`echo $aerm3date | cut -c9-10`

prefix_obs=${cdump}.t${cyc}z
prefix_prep=$prefix_obs
prefix_tbc=gdas1.t${gcyc}z
prefix_inp=${cdump}.t${gcyc}z
prefix_out=${cdump}.t${cyc}z
suffix_obs=${cdump}.${CDATE}
suffix_bias=${cdump}.${gdate}
suffix=tm00.bufr_d
sfcsuffix=nemsio
atmsuffix=nemsio
dtfsuffix="nc4"
DIAG_SUFFIX=".nc4"

aer03=${aerpath}/$m3yy/$m3mm/MERRA2_AER3D_FV3L64.${m3yy}${m3mm}${m3dd}${m3hh}.nc
aer06=${aerpath}/$ayy/$amm/MERRA2_AER3D_FV3L64.${ayy}${amm}${add}${cyc}.nc
aer09=${aerpath}/$p3yy/$p3mm/MERRA2_AER3D_FV3L64.${p3yy}${p3mm}${p3dd}${p3hh}.nc

#datobs=${homedir}/INPUTS/OBSDATA/${cdump}.$PDY/$cyc
datobs=${obsarch}/${cdump}.$PDY/$cyc
datges=${gesarch}/${cdump}.$gPDY/$gcyc

# Set up $tmpdir
rm -rf $tmpdir
mkdir -p $tmpdir
cd $tmpdir
rm -rf core*


# CO2 namelist and file decisions
ICO2=${ICO2:-0}
if [ $ICO2 -gt 0 ] ; then
        # Copy co2 files to $tmpdir
        co2dir=${CO2DIR:-$fixgsi}
        yyyy=$(echo ${CDATE}|cut -c1-4)
        rm ./global_co2_data.txt
        co2=$co2dir/global_co2.gcmscl_$yyyy.txt
        while [ ! -s $co2 ] ; do
                ((yyyy-=1))
                co2=$co2dir/global_co2.gcmscl_$yyyy.txt
        done
        if [ -s $co2 ] ; then
                $ncp $co2 ./global_co2_data.txt
        fi
        if [ ! -s ./global_co2_data.txt ] ; then
                echo "\./global_co2_data.txt" not created
                exit 1
   fi
fi
#CH4 file decision
ICH4=${ICH4:-0}
if [ $ICH4 -gt 0 ] ; then
#        # Copy ch4 files to $tmpdir
        ch4dir=${CH4DIR:-$fixgsi}
        yyyy=$(echo ${CDATE}|cut -c1-4)
        rm ./ch4globaldata.txt
        ch4=$ch4dir/global_ch4_esrlctm_$yyyy.txt
        while [ ! -s $ch4 ] ; do
                ((yyyy-=1))
                ch4=$ch4dir/global_ch4_esrlctm_$yyyy.txt
        done
        if [ -s $ch4 ] ; then
                $ncp $ch4 ./ch4globaldata.txt
        fi
        if [ ! -s ./ch4globaldata.txt ] ; then
                echo "\./ch4globaldata.txt" not created
                exit 1
   fi
fi
IN2O=${IN2O:-0}
if [ $IN2O -gt 0 ] ; then
#        # Copy ch4 files to $tmpdir
        n2odir=${N2ODIR:-$fixgsi}
        yyyy=$(echo ${CDATE}|cut -c1-4)
        rm ./n2oglobaldata.txt
        n2o=$n2odir/global_n2o_esrlctm_$yyyy.txt
        while [ ! -s $n2o ] ; do
                ((yyyy-=1))
                n2o=$n2odir/global_n2o_esrlctm_$yyyy.txt
        done
        if [ -s $n2o ] ; then
                $ncp $n2o ./n2oglobaldata.txt
        fi
        if [ ! -s ./n2oglobaldata.txt ] ; then
                echo "\./n2oglobaldata.txt" not created
                exit 1
   fi
fi
ICO=${ICO:-0}
if [ $ICO -gt 0 ] ; then
#        # Copy CO files to $tmpdir
        codir=${CODIR:-$fixgsi}
        yyyy=$(echo ${CDATE}|cut -c1-4)
        rm ./coglobaldata.txt
        co=$codir/global_co_esrlctm_$yyyy.txt
        while [ ! -s $co ] ; do
                ((yyyy-=1))
                co=$codir/global_co_esrlctm_$yyyy.txt
        done
        if [ -s $co ] ; then
                $ncp $co ./coglobaldata.txt
        fi
        if [ ! -s ./coglobaldata.txt ] ; then
                echo "\./coglobaldata.txt" not created
                exit 1
   fi
fi

# Make gsi namelist

cat << EOF > gsiparm.anl

 &SETUP
   miter=${nummiter},niter(1)=100,niter(2)=50,
   niter_no_qc(1)=2,niter_no_qc(2)=0,
   write_diag(1)=.true.,write_diag(2)=.false.,write_diag(3)=.true.,
   qoption=2,
   gencode=82,factqmin=0.5,factqmax=0.005,deltim=1200,
   iguess=-1,
   tzr_qc=1,
   oneobtest=.false.,retrieval=.false.,l_foto=.false.,
   use_pbl=.false.,use_compress=.true.,nsig_ext=12,gpstop=50.,
   use_gfs_nemsio=.true.,sfcnst_comb=.true.,
   use_readin_anl_sfcmask=.false.,
   lrun_subdirs=.true.,
   crtm_coeffs_path='./crtm_coeffs/',
   newpc4pred=.true.,adp_anglebc=.true.,angord=4,passive_bc=.true.,use_edges=.false.,
   diag_precon=.true.,step_start=1.e-3,emiss_bc=.true.,
   thin4d=.true.,cwoption=3,imp_physics=11,lupp=.true.,
   binary_diag=.false.,netcdf_diag=.true.,
   l4densvar=.false.,ens_nstarthr=3,nhr_obsbin=3,nhr_assimilation=6,lwrite4danl=.false.,
   use_fv3_aero=.true.,
   verbose=${VERBOSE},
   lwrite_peakwt=.true.,
   lread_obs_save=${if_read_obs_save},lread_obs_skip=${if_read_obs_skip},
 /
 &GRIDOPTS
   JCAP_B=$JCAP,JCAP=$JCAP_A,NLAT=$NLAT_A,NLON=$NLON_A,nsig=$LEVS,
   regional=.false.,nlayers(63)=3,nlayers(64)=6,

 /
 &BKGERR
   vs=0.7,
   hzscl=1.7,0.8,0.5,
   hswgt=0.45,0.3,0.25,
   bw=0.0,norsp=4,
   bkgv_flowdep=.true.,bkgv_rewgtfct=1.5,
   bkgv_write=.false.,
   cwcoveqqcov=.false.,

 /
 &ANBKGERR
   anisotropic=.false.,

 /
 &JCOPTS
   ljcdfi=.false.,alphajc=0.0,ljcpdry=.true.,bamp_jcpdry=2.5e7,ljc4tlevs=.true.,

 /
 &STRONGOPTS
   tlnmc_option=2,nstrong=1,nvmodes_keep=8,period_max=6.,period_width=1.5,
   baldiag_full=.false.,baldiag_inc=.false.,

 /
 &OBSQC
   dfact=0.75,dfact1=3.0,noiqc=.true.,oberrflg=.false.,c_varqc=0.04,
   use_poq7=.true.,qc_noirjaco3_pole=.true.,vqc=.true.,
   aircraft_t_bc=.true.,biaspredt=1000.0,upd_aircraft=.true.,cleanup_tail=.true.

 /
 &OBS_INPUT
   dmesh(1)=145.0,dmesh(2)=150.0,dmesh(3)=100.0,time_window_max=3.0,

 /
OBS_INPUT::
!  dfile          dtype       dplat       dsis                dval    dthin dsfcalc
!   prepbufr       ps          null        ps                  0.0     0     0
!   prepbufr       t           null        t                   0.0     0     0
!   prepbufr_profl t           null        t                   0.0     0     0
!   prepbufr       q           null        q                   0.0     0     0
!   prepbufr_profl q           null        q                   0.0     0     0
!   prepbufr       pw          null        pw                  0.0     0     0
!   prepbufr       uv          null        uv                  0.0     0     0
!   prepbufr_profl uv          null        uv                  0.0     0     0
!   satwndbufr     uv          null        uv                  0.0     0     0
!   prepbufr       spd         null        spd                 0.0     0     0
!   prepbufr       dw          null        dw                  0.0     0     0
!   radarbufr      rw          null        rw                  0.0     0     0
!   nsstbufr       sst         nsst        sst                 0.0     0     0
!   gpsrobufr      gps_bnd     null        gps                 0.0     0     0
!   ssmirrbufr     pcp_ssmi    dmsp        pcp_ssmi            0.0    -1     0
!   tmirrbufr      pcp_tmi     trmm        pcp_tmi             0.0    -1     0
!   sbuvbufr       sbuv2       n16         sbuv8_n16           0.0     0     0
!   sbuvbufr       sbuv2       n17         sbuv8_n17           0.0     0     0
!   sbuvbufr       sbuv2       n18         sbuv8_n18           0.0     0     0
!   hirs3bufr      hirs3       n17         hirs3_n17           0.0     1     0
!   hirs4bufr      hirs4       metop-a     hirs4_metop-a       0.0     1     1
!   gimgrbufr      goes_img    g11         imgr_g11            0.0     1     0
!   gimgrbufr      goes_img    g12         imgr_g12            0.0     1     0
   airsbufr       airs        aqua        airs_aqua           0.0     1     1
!   amsuabufr      amsua       n15         amsua_n15           0.0     1     1
!   amsuabufr      amsua       n18         amsua_n18           0.0     1     1
!   amsuabufr      amsua       metop-a     amsua_metop-a       0.0     1     1
!   airsbufr       amsua       aqua        amsua_aqua          0.0     1     1
!   amsubbufr      amsub       n17         amsub_n17           0.0     1     1
!   mhsbufr        mhs         n18         mhs_n18             0.0     1     1
!   mhsbufr        mhs         metop-a     mhs_metop-a         0.0     1     1
!   ssmitbufr      ssmi        f15         ssmi_f15            0.0     1     0
!   amsrebufr      amsre_low   aqua        amsre_aqua          0.0     1     0
!   amsrebufr      amsre_mid   aqua        amsre_aqua          0.0     1     0
!   amsrebufr      amsre_hig   aqua        amsre_aqua          0.0     1     0
!   ssmisbufr      ssmis       f16         ssmis_f16           0.0     1     0
!   ssmisbufr      ssmis       f17         ssmis_f17           0.0     1     0
!   ssmisbufr      ssmis       f18         ssmis_f18           0.0     1     0
!   ssmisbufr      ssmis       f19         ssmis_f19           0.0     1     0
!   gsnd1bufr      sndrd1      g12         sndrD1_g12          0.0     1     0
!   gsnd1bufr      sndrd2      g12         sndrD2_g12          0.0     1     0
!   gsnd1bufr      sndrd3      g12         sndrD3_g12          0.0     1     0
!   gsnd1bufr      sndrd4      g12         sndrD4_g12          0.0     1     0
!   gsnd1bufr      sndrd1      g11         sndrD1_g11          0.0     1     0
!   gsnd1bufr      sndrd2      g11         sndrD2_g11          0.0     1     0
!   gsnd1bufr      sndrd3      g11         sndrD3_g11          0.0     1     0
!   gsnd1bufr      sndrd4      g11         sndrD4_g11          0.0     1     0
!   gsnd1bufr      sndrd1      g13         sndrD1_g13          0.0     1     0
!   gsnd1bufr      sndrd2      g13         sndrD2_g13          0.0     1     0
!   gsnd1bufr      sndrd3      g13         sndrD3_g13          0.0     1     0
!   gsnd1bufr      sndrd4      g13         sndrD4_g13          0.0     1     0
   iasibufr       iasi        metop-a     iasi_metop-a        0.0     1     1
!   gomebufr       gome        metop-a     gome_metop-a        0.0     2     0
!   omibufr        omi         aura        omi_aura            0.0     2     0
!   sbuvbufr       sbuv2       n19         sbuv8_n19           0.0     0     0
!   hirs4bufr      hirs4       n19         hirs4_n19           0.0     1     1
!   amsuabufr      amsua       n19         amsua_n19           0.0     1     1
!   mhsbufr        mhs         n19         mhs_n19             0.0     1     1
!   tcvitl         tcp         null        tcp                 0.0     0     0
!   seviribufr     seviri      m08         seviri_m08          0.0     1     0
!   seviribufr     seviri      m09         seviri_m09          0.0     1     0
!   seviribufr     seviri      m10         seviri_m10          0.0     1     0
!   hirs4bufr      hirs4       metop-b     hirs4_metop-b       0.0     1     1
!   amsuabufr      amsua       metop-b     amsua_metop-b       0.0     1     1
!   mhsbufr        mhs         metop-b     mhs_metop-b         0.0     1     1
   iasibufr       iasi        metop-b     iasi_metop-b        0.0     1     1
!   gomebufr       gome        metop-b     gome_metop-b        0.0     2     0
!   atmsbufr       atms        npp         atms_npp            0.0     1     1
!   atmsbufr       atms        n20         atms_n20            0.0     1     1
   crisbufr       cris        npp         cris_npp            0.0     1     0
   crisfsbufr     cris-fsr    npp         cris-fsr_npp        0.0     1     0
   crisfsbufr     cris-fsr    n20         cris-fsr_n20        0.0     1     0
!   gsnd1bufr      sndrd1      g14         sndrD1_g14          0.0     1     0
!   gsnd1bufr      sndrd2      g14         sndrD2_g14          0.0     1     0
!   gsnd1bufr      sndrd3      g14         sndrD3_g14          0.0     1     0
!   gsnd1bufr      sndrd4      g14         sndrD4_g14          0.0     1     0
!   gsnd1bufr      sndrd1      g15         sndrD1_g15          0.0     1     0
!   gsnd1bufr      sndrd2      g15         sndrD2_g15          0.0     1     0
!   gsnd1bufr      sndrd3      g15         sndrD3_g15          0.0     1     0
!   gsnd1bufr      sndrd4      g15         sndrD4_g15          0.0     1     0
!   oscatbufr      uv          null        uv                  0.0     0     0
!   mlsbufr        mls30       aura        mls30_aura          0.0     0     0
!   avhambufr      avhrr       metop-a     avhrr3_metop-a      0.0     1     0
!   avhpmbufr      avhrr       n18         avhrr3_n18          0.0     1     0
!   amsr2bufr      amsr2       gcom-w1     amsr2_gcom-w1       0.0     3     0
!   gmibufr        gmi         gpm         gmi_gpm             0.0     3     0
!   saphirbufr     saphir      meghat      saphir_meghat       0.0     3     0
!   ahibufr        ahi         himawari8   ahi_himawari8       0.0     3     0
!   rapidscatbufr  uv          null        uv                  0.0     0     0
::
  &SUPEROB_RADAR

 /
 &LAG_DATA

 /
 &HYBRID_ENSEMBLE
   l_hyb_ens=.false.,n_ens=20,beta_s0=0.125,readin_beta=.false.,s_ens_h=800,s_ens_v=-0.8,generate_ens=.false.,uv_hyb_ens=.true.,jcap_ens=62,
   nlat_ens=96,nlon_ens=192,ANISO_A_EN=.false.,jcap_ens_test=62,oz_univ_static=.false.,readin_localization=.true.,ensemble_path='./ensemble_data/',
   ens_fast_read=.true.,write_ens_sprd=.false.,

 /
 &RAPIDREFRESH_CLDSURF
   dfi_radar_latent_heat_time_period=30.0,
 /
 &CHEM
   lread_ext_aerosol=${READEXTAER},lmerra2aer=${MERRA2AER}
 /
 &SINGLEOB_TEST
   ${SINGLEOB_TEST}
 /
 &NST
   nst_gsi=3,nstinfo=4,fac_dtl=1,fac_tsl=1,zsea1=0,zsea2=0,

 /
EOF

# Set fixed files
#   berror   = forecast model background error statistics
#   specoef  = CRTM spectral coefficients
#   trncoef  = CRTM transmittance coefficients
#   emiscoef = CRTM coefficients for IR sea surface emissivity model
#   aerocoef = CRTM coefficients for aerosol effects
#   cldcoef  = CRTM coefficients for cloud effects
#   satinfo  = text file with information about assimilation of brightness temperatures
#   cloudyinfo  = text file with information about assimilation of cloudy radiance
#   satangl  = angle dependent bias correction file (fixed in time)
#   pcpinfo  = text file with information about assimilation of prepcipitation rates
#   ozinfo   = text file with information about assimilation of ozone data
#   errtable = text file with obs error for conventional data (optional)
#   convinfo = text file with information about assimilation of conventional data
#   bufrtable= text file ONLY needed for single obs test (oneobstest=.true.)
#   bftab_sst= bufr table for sst ONLY needed for sst retrieval (retrieval=.true.)

berror=$fixgsi/Big_Endian/global_berror.l${LEVS}y${NLAT_A}.f77

emiscoef_IRwater=$fixcrtm/Nalli.IRwater.EmisCoeff.bin
emiscoef_IRice=$fixcrtm/NPOESS.IRice.EmisCoeff.bin
emiscoef_IRland=$fixcrtm/NPOESS.IRland.EmisCoeff.bin
emiscoef_IRsnow=$fixcrtm/NPOESS.IRsnow.EmisCoeff.bin
emiscoef_VISice=$fixcrtm/NPOESS.VISice.EmisCoeff.bin
emiscoef_VISland=$fixcrtm/NPOESS.VISland.EmisCoeff.bin
emiscoef_VISsnow=$fixcrtm/NPOESS.VISsnow.EmisCoeff.bin
emiscoef_VISwater=$fixcrtm/NPOESS.VISwater.EmisCoeff.bin
emiscoef_MWwater=$fixcrtm/FASTEM6.MWwater.EmisCoeff.bin
aercoef=$fixcrtm/AerosolCoeff.bin
cldcoef=$fixcrtm/CloudCoeff.bin
satangl=$fixgsi/global_satangbias.txt
scaninfo=$fixgsi/global_scaninfo.txt
cloudyinfo=$fixgsi/cloudy_radiance_info.txt
convinfo=$fixgsi/global_convinfo_reg_test.txt
aeroinfo=$fixgsi/aeroinfo_fv3aerorad
ozinfo=$fixgsi/global_ozinfo.txt
pcpinfo=$fixgsi/global_pcpinfo.txt
hybens_info=$fixgsi/global_hybens_info.l64.txt
errtable=$fixgsi/prepobs_errtable.global
atmsbeaminfo=$fixgsi/atms_beamwidth.txt

# Copy executable and fixed files to $tmpdir
   $ncp $gsiexec  ./gsi.x

mkdir ./crtm_coeffs
$ncp $berror   ./berror_stats
$ncp $emiscoef_IRwater ./crtm_coeffs/Nalli.IRwater.EmisCoeff.bin
$ncp $emiscoef_IRice ./crtm_coeffs/NPOESS.IRice.EmisCoeff.bin
$ncp $emiscoef_IRsnow ./crtm_coeffs/NPOESS.IRsnow.EmisCoeff.bin
$ncp $emiscoef_IRland ./crtm_coeffs/NPOESS.IRland.EmisCoeff.bin
$ncp $emiscoef_VISice ./crtm_coeffs/NPOESS.VISice.EmisCoeff.bin
$ncp $emiscoef_VISland ./crtm_coeffs/NPOESS.VISland.EmisCoeff.bin
$ncp $emiscoef_VISsnow ./crtm_coeffs/NPOESS.VISsnow.EmisCoeff.bin
$ncp $emiscoef_VISwater ./crtm_coeffs/NPOESS.VISwater.EmisCoeff.bin
$ncp $emiscoef_MWwater ./crtm_coeffs/FASTEM6.MWwater.EmisCoeff.bin
$ncp $aercoef  ./crtm_coeffs/AerosolCoeff.bin
$ncp $cldcoef  ./crtm_coeffs/CloudCoeff.bin
$ncp $satangl  ./satbias_angle
$ncp $scaninfo ./scaninfo
$ncp $satinfo  ./satinfo
$ncp $cloudyinfo  ./cloudy_radiance_info.txt
$ncp $pcpinfo  ./pcpinfo
$ncp $ozinfo   ./ozinfo
$ncp $convinfo ./convinfo
$ncp $errtable ./errtable
$ncp $anavinfo ./anavinfo
$ncp $aeroinfo ./aeroinfo
$ncp $hybens_info ./hybens_info
$ncp $atmsbeaminfo ./atms_beamwidth.txt

# Copy CRTM coefficient files
for file in `awk '{if($1!~"!"){print $1}}' ./satinfo | sort | uniq` ;do
    $ncp $fixcrtm/${file}.SpcCoeff.bin ./crtm_coeffs
    $ncp $fixcrtm/${file}.TauCoeff.bin ./crtm_coeffs
done

# Copy observational data to $tmpdir
$ncpl $datobs/${prefix_obs}.prepbufr                ./prepbufr
$ncpl $datobs/${prefix_obs}.prepbufr.acft_profiles  ./prepbufr_profl
if [ $machine == 'hera' ]; then
   $ncpl $datobs/${prefix_obs}.nsstbufr                ./nsstbufr
elif [ $machine == 's4' ]; then
   cat $datobs/${prefix_obs}.sfcshp.${suffix} \
       $datobs/${prefix_obs}.tesac.$suffix \
       $datobs/${prefix_obs}.bathy.$suffix \
       $datobs/${prefix_obs}.trkob.$suffix > ./nsstbufr
fi
$ncpl $datobs/${prefix_obs}.gpsro.${suffix}         ./gpsrobufr
$ncpl $datobs/${prefix_obs}.satwnd.${suffix}        ./satwndbufr
$ncpl $datobs/${prefix_obs}.spssmi.${suffix}        ./ssmirrbufr
$ncpl $datobs/${prefix_obs}.sptrmm.${suffix}        ./tmirrbufr
$ncpl $datobs/${prefix_obs}.osbuv8.${suffix}        ./sbuvbufr
$ncpl $datobs/${prefix_obs}.goesfv.${suffix}        ./gsnd1bufr
$ncpl $datobs/${prefix_obs}.1bamua.${suffix}        ./amsuabufr
$ncpl $datobs/${prefix_obs}.1bamub.${suffix}        ./amsubbufr
$ncpl $datobs/${prefix_obs}.1bhrs2.${suffix}        ./hirs2bufr
$ncpl $datobs/${prefix_obs}.1bhrs3.${suffix}        ./hirs3bufr
$ncpl $datobs/${prefix_obs}.1bhrs4.${suffix}        ./hirs4bufr
$ncpl $datobs/${prefix_obs}.1bmhs.${suffix}         ./mhsbufr
$ncpl $datobs/${prefix_obs}.1bmsu.${suffix}         ./msubufr
$ncpl $datobs/${prefix_obs}.airsev.${suffix}        ./airsbufr
$ncpl $datobs/${prefix_obs}.sevcsr.${suffix}        ./seviribufr
$ncpl $datobs/${prefix_obs}.mtiasi.${suffix}        ./iasibufr
$ncpl $datobs/${prefix_obs}.ssmit.${suffix}         ./ssmitbufr
$ncpl $datobs/${prefix_obs}.ssmisu.${suffix}        ./ssmisbufr
$ncpl $datobs/${prefix_obs}.gome.${suffix}          ./gomebufr
$ncpl $datobs/${prefix_obs}.omi.${suffix}           ./omibufr
$ncpl $datobs/${prefix_obs}.mls.${suffix}           ./mlsbufr
$ncpl $datobs/${prefix_obs}.eshrs3.${suffix}        ./hirs3bufrears
$ncpl $datobs/${prefix_obs}.esamua.${suffix}        ./amsuabufrears
$ncpl $datobs/${prefix_obs}.esamub.${suffix}        ./amsubbufrears
$ncpl $datobs/${prefix_obs}.atms.${suffix}          ./atmsbufr
$ncpl $datobs/${prefix_obs}.cris.${suffix}          ./crisbufr
$ncpl $datobs/${prefix_obs}.crisf4.${suffix}        ./crisfsbufr
$ncpl $datobs/${prefix_obs}.syndata.tcvitals.tm00   ./tcvitl
$ncpl $datobs/${prefix_obs}.avcsam.${suffix}        ./avhambufr
$ncpl $datobs/${prefix_obs}.avcspm.${suffix}        ./avhpmbufr
$ncpl $datobs/${prefix_obs}.saphir.${suffix}        ./saphirbufr
$ncpl $datobs/${prefix_obs}.gmi1cr.${suffix}        ./gmibufr
if [ "$debug" = ".false." ]; then
    $ncpl $datobs/${prefix_obs}.esiasi.${suffix}        ./iasibufrears
fi
$ncpl $datobs/${prefix_obs}.hrs3db.${suffix}        ./hirs3bufr_db
$ncpl $datobs/${prefix_obs}.amuadb.${suffix}        ./amsuabufr_db
$ncpl $datobs/${prefix_obs}.amubdb.${suffix}        ./amsubbufr_db
$ncpl $datobs/${prefix_obs}.iasidb.${suffix}        ./iasibufr_db
$ncpl $datobs/${prefix_obs}.crisdb.${suffix}        ./crisbufr_db
$ncpl $datobs/${prefix_obs}.atmsdb.${suffix}        ./atmsbufr_db
$ncpl $datobs/${prefix_obs}.escris.${suffix}        ./crisbufrears
$ncpl $datobs/${prefix_obs}.esatms.${suffix}        ./atmsbufrears


# Copy bias correction, atmospheric and surface files
$ncpl $datges/${cdump}.t${gcyc}z.abias           ./satbias_in
$ncpl $datges/${cdump}.t${gcyc}z.abias_pc        ./satbias_pc
$ncpl $datges/${cdump}.t${gcyc}z.abias_air       ./aircftbias_in
#$ncpl $datges/gfs.t18z.radstat         ./radstat.gdas

$ncpl $aer03  ./aerf03
$ncpl $aer06  ./aerf06
$ncpl $aer09  ./aerf09

$ncpl $datges/${prefix_inp}.sfcf003.${sfcsuffix}       ./sfcf03
$ncpl $datges/${prefix_inp}.sfcf006.${sfcsuffix}       ./sfcf06
$ncpl $datges/${prefix_inp}.sfcf009.${sfcsuffix}       ./sfcf09

$ncpl $datges/${prefix_inp}.atmf003.${atmsuffix}       ./sigf03
$ncpl $datges/${prefix_inp}.atmf006.${atmsuffix}       ./sigf06
$ncpl $datges/${prefix_inp}.atmf009.${atmsuffix}       ./sigf09

# Run GSI
echo "run gsi now"
eval "$APRUN ./gsi.x 2>&1 | tee stdout"
#exit
rc=$PIPESTATUS
if [ $rc -eq 0 ]; then
   mv siganl              $outdir/${prefix_out}.atmanl.${atmsuffix}
   mv dtfanl              $outdir/${prefix_out}.dtfanl.${dtfsuffix}
   mv stdout              $outdir/gsi.stdout.${CDATE}
else
   echo "GSI failed!!"
   exit 1
fi


# Loop over first and last outer loops to generate innovation
# diagnostic files for indicated observation types (groups)
#
# NOTE:  Since we set miter=2 in GSI namelist SETUP, outer
#        loop 03 will contain innovations with respect to 
#        the analysis.  Creation of o-a innovation files
#        is triggered by write_diag(3)=.true.  The setting
#        write_diag(1)=.true. turns on creation of o-g
#        innovation files.
#

echo "Time before diagnostic loop is `date` "

ntype=3

diagtype[0]="conv conv_gps conv_ps conv_pw conv_q conv_sst conv_t conv_tcp conv_uv conv_spd"
diagtype[1]="pcp_ssmi_dmsp pcp_tmi_trmm"
diagtype[2]="sbuv2_n16 sbuv2_n17 sbuv2_n18 sbuv2_n19 gome_metop-a gome_metop-b omi_aura mls30_aura ompsnp_npp ompstc8_npp gome_metop-c"
diagtype[3]="hirs2_n14 msu_n14 sndr_g08 sndr_g11 sndr_g12 sndr_g13 sndr_g08_prep sndr_g11_prep sndr_g12_prep sndr_g13_prep sndrd1_g11 sndrd2_g11 sndrd3_g11 sndrd4_g11 sndrd1_g12 sndrd2_g12 sndrd3_g12 sndrd4_g12 sndrd1_g13 sndrd2_g13 sndrd3_g13 sndrd4_g13 sndrd1_g14 sndrd2_g14 sndrd3_g14 sndrd4_g14 sndrd1_g15 sndrd2_g15 sndrd3_g15 sndrd4_g15 hirs3_n15 hirs3_n16 hirs3_n17 amsua_n15 amsua_n16 amsua_n17 amsub_n15 amsub_n16 amsub_n17 hsb_aqua airs_aqua amsua_aqua imgr_g08 imgr_g11 imgr_g12 imgr_g14 imgr_g15 ssmi_f13 ssmi_f15 hirs4_n18 hirs4_metop-a amsua_n18 amsua_metop-a mhs_n18 mhs_metop-a amsre_low_aqua amsre_mid_aqua amsre_hig_aqua ssmis_f16 ssmis_f17 ssmis_f18 ssmis_f19 ssmis_f20 iasi_metop-a hirs4_n19 amsua_n19 mhs_n19 seviri_m08 seviri_m09 seviri_m10 seviri_m11 cris_npp cris-fsr_npp cris-fsr_n20 atms_npp atms_n20 hirs4_metop-b amsua_metop-b mhs_metop-b iasi_metop-b avhrr_metop-b avhrr_n18 avhrr_n19 avhrr_metop-a amsr2_gcom-w1 gmi_gpm saphir_meghat ahi_himawari8 abi_g16 abi_g17 amsua_metop-c mhs_metop-c iasi_metop-c avhrr_metop-c"

numfile[0]=0
numfile[1]=0
numfile[2]=0
numfile[3]=0

# Set diagnostic file prefix based on lrun_subdirs variable
prefix=" dir.*/"

cd $tmpdir

loops="01 03"
for loop in $loops; do
   case $loop in
      01) string=ges;;
      03) string=anl;;
       *) string=$loop;;
   esac
   echo $(date) START loop $string >&2
   n=-1
   while [ $((n+=1)) -le $ntype ] ;do
      for type in $(echo ${diagtype[n]}); do
         count=$(ls ${prefix}${type}_${loop}* 2>/dev/null | wc -l)
         if [ $count -gt 1 ]; then
            if [ $binary_diag = ".true." ]; then
               cat ${prefix}${type}_${loop}* > diag_${type}_${string}.${CDATE}${DIAG_SUFFIX}
            else
               $CATEXEC -o diag_${type}_${string}.${CDATE}${DIAG_SUFFIX} ${prefix}${type}_${loop}*
            fi
            numfile[n]=$(expr ${numfile[n]} + 1)
         elif [ $count -eq 1 ]; then
            cat ${prefix}${type}_${loop}* > diag_${type}_${string}.${CDATE}${DIAG_SUFFIX}
            numfile[n]=$(expr ${numfile[n]} + 1)
         fi
         if [ -s $tmpdir/diag_${type}_${string}.${CDATE}${DIAG_SUFFIX} ]; then
            mv $tmpdir/diag_${type}_${string}.${CDATE}${DIAG_SUFFIX} $outdir
         fi
      done
   done
   echo $(date) END loop $string >&2
done

# If requested, clean up $tmpdir
if [[ "$CLEAN" = "YES" ]];then
   if [[ $rc -eq 0 ]];then
      rm -rf $tmpdir
      cd $tmpdir
      cd ../
      rmdir $tmpdir
   fi
fi


# End of script
exit
