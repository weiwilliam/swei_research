#!/bin/sh --login
#SBATCH --output=../logs/idealcrtm.%j
#SBATCH --job-name=swei_idealcrtm
#SBATCH --qos=batch
##SBATCH --qos=debug
#SBATCH --time=1:00:00
#SBATCH --nodes=1 
#SBATCH --account=gsd-fv3-test
#

machine='s4'

if [ $machine == 'hera' ]; then
   HOMEDIR=/scratch1/BMC/gsd-fv3-dev/Shih-wei.Wei/IdealCRTM
   WRKDIR=${HOMEDIR}/wrk
   USHDIR=${HOMEDIR}/ush
   OUTDIR=${HOMEDIR}/output
   EXEC=${HOMEDIR}/bin/crtm_ideal_FDM
   # Module load
   . /apps/lmod/lmod/init/sh
   module purge
   module load intel
   module load netcdf
elif [ $machine == 's4' ]; then
   SCRPDIR=/home/swei/research/IdealCRTM
   USHDIR=${SCRPDIR}/ush
   HOMEDIR=/data/users/swei/Experiments/IdealCRTM
   WRKDIR=${HOMEDIR}/wrk
   OUTDIR=${HOMEDIR}/output
   EXEC=${HOMEDIR}/bin/crtm_ideal_FDM
   module purge
   module load license_intel/S4
   module load intel
   module load hdf hdf5
   module load netcdf4
fi
module list
INNML=${USHDIR}/idealcrtm.nml.IN
KMODFILE='kmod_inspection.nc'
AODKFILE='aodk_inspection.nc'

# Exp on/off
BSL='N'  # No aerosols experiment
EXP1='Y' # Concentration
EXP2='N' # Thickness (only for genmethod 1)
EXP3='Y' # Altitude
EXP4='Y' # Sfc_Peak_Ratio (SPR), available for genmethod 2 only
EXP5='Y' # Bins Partition
EXP6='N' # Surface type
EXP7='N' # Surface Emissivity

# Default setting for each factor
naers=5
useremi='.false.'
aername="gocart_dust"
case $aername in
'gocart_dust')
   binlist='DU001,DU002,DU003,DU004,DU005'
   binspar='0.1,0.4,0.3,0.15,0.05'
   totconc='1.967e-03'
   binsparlst="0.1,0.4,0.3,0.15,0.05 0.15,0.45,0.35,0.05,0."
   conclist="2.45875e-04 4.9175e-04 9.835e-04 1.967e-03 3.934e-03 5.901e-03"
   ;;
'gocart_carbon')
   binlist='ocphobic,ocphilic,bcphobic,bcphilic'
   binspar='0.15,0.70,0.05,0.1'
   totconc='1.818e-04'
   binsparlst="0.15,0.7,0.05,0.1 0.05,0.9,0.,0.05"
   conclist='2.2725e-05 4.5450e-05 9.0900e-05 1.8180e-04 3.6360e-04 5.4540e-04'
   ;;
'gocart_seas')
   binlist='SS001,SS002,SS003,SS004'
   binspar='0.05,0.25,0.6,0.1'
   totconc='7.176e-04'
   conclist='8.9700e-05 1.7940e-04 3.5880e-04 7.1760e-04 1.4352e-03 2.1528e-03'
   ;;
'gocart_sulf')
   binlist='sulfate'
   binspar='1.'
   totconc='1.2527e-04'
   conclist='1.5659e-05 3.1318e-05 6.2635e-05 1.2527e-04 2.5054e-04 3.7581e-04'
   ;;
esac
genmethod='2'
landcover='0.'
landtype='10'
lai='0.17'

#genmethod=1
thickns='1'
atlevel='76'
thicklist="1 3 5 7 9"

# genmethod=2
modelv='76'
sfc_mode_ratio='0.5'
bmoderatio='0.8'

daytime='.false.'

atlvllist="76 83 88 92"
smrlist="0.8 0.5 0.2 0"
sfctypelst="water desert"

rm -rf $WRKDIR/*
if [ ! -s $WRKDIR/coefficients ]; then
   mkdir -p $WRKDIR/coefficients
   cd $WRKDIR/coefficients
   sh $USHDIR/link_crtm_coeff.sh $machine
fi
cd $WRKDIR

case $BSL in
Y|y)
echo 'no aerosol baseline run over different surface'
for sfc in $sfctypelst
do
  case $sfc in
  'desert')
     explandcover='1.0' ; expsfctype='3'  ; explai='0.0' ;;
  'tundra')
     explandcover='1.0' ; expsfctype='10' ; explai='0.17';;
  'water')
     explandcover='0.0' ; expsfctype='10' ; explai='0.17';;
  *)
    echo 'not available surface type test' ;  exit ;;
  esac
  echo $sfc
  [[ -s idealcrtm.nml ]]&&rm idealcrtm.nml
  cat ${INNML} | sed s:_NAERS_:0:g \
               | sed s:_USEREMI_:${useremi}:g \
               | sed s:_AERNAME_:${aername}:g \
               | sed s:_BINSLST_:${binlist}:g \
               | sed s:_BINSPAR_:${binspar}:g \
               | sed s:_GENMETHOD_:${genmethod}:g \
               | sed s:_TOTCONC_:${totconc}:g \
               | sed s:_THICKNS_:${thickns}:g \
               | sed s:_ATLEVEL_:${atlevel}:g \
               | sed s:_MODELV_:${modelv}:g \
               | sed s:_BMODERATIO_:${bmoderatio}:g \
               | sed s:_SFCMODRATIO_:${sfc_mode_ratio}:g \
               | sed s:_LANDCOV_:${explandcover}:g \
               | sed s:_LANDTYPE_:${expsfctype}:g \
               | sed s:_LAI_:${explai}:g \
               | sed s:\'::g \
               > idealcrtm.nml
  if [ -s idealcrtm.nml ] ; then
     $EXEC
     rc=$?
     if [ $rc -eq 0 ]; then
        mkdir -p $OUTDIR/NoAer
        mv $KMODFILE $OUTDIR/NoAer/noaer_${sfc}_kmod.nc
        mv $AODKFILE $OUTDIR/NoAer/noaer_${sfc}_aodk.nc
     elif [ $rc -ne 0 ]; then
        exit 1
     fi
  fi
done
;;
*)
echo 'Skip clear-sky test' ;;
esac


case $EXP1 in
Y|y)
# Total Concentration Sensitivity
for conc in $conclist
do
  echo $conc
  [[ -s idealcrtm.nml ]]&&rm idealcrtm.nml
  cat ${INNML} | sed s:_NAERS_:${naers}:g \
               | sed s:_USEREMI_:${useremi}:g \
               | sed s:_AERNAME_:${aername}:g \
               | sed s:_BINSLST_:${binlist}:g \
               | sed s:_BINSPAR_:${binspar}:g \
               | sed s:_GENMETHOD_:${genmethod}:g \
               | sed s:_TOTCONC_:${conc}:g \
               | sed s:_THICKNS_:${thickns}:g \
               | sed s:_ATLEVEL_:${atlevel}:g \
               | sed s:_MODELV_:${modelv}:g \
               | sed s:_BMODERATIO_:${bmoderatio}:g \
               | sed s:_SFCMODRATIO_:${sfc_mode_ratio}:g \
               | sed s:_LANDCOV_:${landcover}:g \
               | sed s:_LANDTYPE_:${landtype}:g \
               | sed s:_LAI_:${lai}:g \
               | sed s:\'::g \
               > idealcrtm.nml
  if [ -s idealcrtm.nml ] ; then
     $EXEC
     rc=$?
     if [ $rc -eq 0 ]; then
        mkdir -p $OUTDIR/Concentration
        mv $KMODFILE $OUTDIR/Concentration/${conc}_kmod.nc
        mv $AODKFILE $OUTDIR/Concentration/${conc}_aodk.nc
     elif [ $rc -ne 0 ]; then
        exit 1
     fi
  fi
done
;;
*)
echo 'Skip concentration test' ;;
esac

case $EXP2 in
Y|y)
# Thickness Sensitivity
for thick in $thicklist
do
  echo $thick
  [[ -s idealcrtm.nml ]]&&rm idealcrtm.nml
  cat ${INNML} | sed s:_NAERS_:${naers}:g \
               | sed s:_USEREMI_:${useremi}:g \
               | sed s:_AERNAME_:${aername}:g \
               | sed s:_BINSLST_:${binlist}:g \
               | sed s:_BINSPAR_:${binspar}:g \
               | sed s:_GENMETHOD_:${genmethod}:g \
               | sed s:_TOTCONC_:${totconc}:g \
               | sed s:_THICKNS_:${thick}:g \
               | sed s:_ATLEVEL_:${atlevel}:g \
               | sed s:_MODELV_:${modelv}:g \
               | sed s:_BMODERATIO_:${bmoderatio}:g \
               | sed s:_SFCMODRATIO_:${sfc_mode_ratio}:g \
               | sed s:_LANDCOV_:${landcover}:g \
               | sed s:_LANDTYPE_:${landtype}:g \
               | sed s:_LAI_:${lai}:g \
               | sed s:\'::g \
               > idealcrtm.nml
  if [ -s idealcrtm.nml ] ; then
     $EXEC
     rc=$?
     if [ $rc -eq 0 ]; then
        mkdir -p $OUTDIR/Thickness
        mv $KMODFILE $OUTDIR/Thickness/L${thick}_at_L${atlevel}_kmod.nc
        mv $AODKFILE $OUTDIR/Thickness/L${thick}_at_L${atlevel}_aodk.nc
     elif [ $rc -ne 0 ]; then
        exit 1
     fi
  fi
done
;;
*)
echo 'Skip thickness test' ;;
esac

case $EXP3 in
Y|y)
# Altitude Sensitivity
for atlvl in $atlvllist
do
  echo $atlvl
  [[ -s idealcrtm.nml ]]&&rm idealcrtm.nml
  cat ${INNML} | sed s:_NAERS_:${naers}:g \
               | sed s:_USEREMI_:${useremi}:g \
               | sed s:_AERNAME_:${aername}:g \
               | sed s:_BINSLST_:${binlist}:g \
               | sed s:_BINSPAR_:${binspar}:g \
               | sed s:_GENMETHOD_:${genmethod}:g \
               | sed s:_TOTCONC_:${totconc}:g \
               | sed s:_THICKNS_:${thickns}:g \
               | sed s:_ATLEVEL_:${atlvl}:g \
               | sed s:_MODELV_:${atlvl}:g \
               | sed s:_BMODERATIO_:${bmoderatio}:g \
               | sed s:_SFCMODRATIO_:${sfc_mode_ratio}:g \
               | sed s:_LANDCOV_:${landcover}:g \
               | sed s:_LANDTYPE_:${landtype}:g \
               | sed s:_LAI_:${lai}:g \
               | sed s:\'::g \
               > idealcrtm.nml
  if [ -s idealcrtm.nml ] ; then
     $EXEC
     rc=$?
     if [ $rc -eq 0 ]; then
        mkdir -p $OUTDIR/Altitude
        if [ $genmethod -eq 1 ]; then
           mv $KMODFILE $OUTDIR/Altitude/L${thickns}_at_L${atlvl}_kmod.nc
           mv $AODKFILE $OUTDIR/Altitude/L${thickns}_at_L${atlvl}_aodk.nc
        elif [ $genmethod -eq 2 ]; then
           mv $KMODFILE $OUTDIR/Altitude/Peak_at_L${atlvl}_kmod.nc
           mv $AODKFILE $OUTDIR/Altitude/Peak_at_L${atlvl}_aodk.nc
        fi 
     elif [ $rc -ne 0 ]; then
        exit 1
     fi
  fi
done
;;
*)
echo 'Skip altitude test' ;;
esac

case $EXP4 in
Y|y)
# Surface Peak Ratio Sensitivity
for smr in $smrlist
do
  echo $smr
  [[ -s idealcrtm.nml ]]&&rm idealcrtm.nml
  cat ${INNML} | sed s:_NAERS_:${naers}:g \
               | sed s:_USEREMI_:${useremi}:g \
               | sed s:_AERNAME_:${aername}:g \
               | sed s:_BINSLST_:${binlist}:g \
               | sed s:_BINSPAR_:${binspar}:g \
               | sed s:_GENMETHOD_:${genmethod}:g \
               | sed s:_TOTCONC_:${totconc}:g \
               | sed s:_THICKNS_:${thickns}:g \
               | sed s:_ATLEVEL_:${atlvl}:g \
               | sed s:_MODELV_:${modelv}:g \
               | sed s:_BMODERATIO_:${bmoderatio}:g \
               | sed s:_SFCMODRATIO_:${smr}:g \
               | sed s:_LANDCOV_:${landcover}:g \
               | sed s:_LANDTYPE_:${landtype}:g \
               | sed s:_LAI_:${lai}:g \
               | sed s:\'::g \
               > idealcrtm.nml
  if [ -s idealcrtm.nml ] ; then
     $EXEC
     rc=$?
     if [ $rc -eq 0 ]; then
        mkdir -p $OUTDIR/Sfc_Peak_Ratio
        mv $KMODFILE $OUTDIR/Sfc_Peak_Ratio/SPR${smr}_kmod.nc
        mv $AODKFILE $OUTDIR/Sfc_Peak_Ratio/SPR${smr}_aodk.nc
     elif [ $rc -ne 0 ]; then
        exit 1
     fi
  fi
done
;;
*)
echo 'Skip altitude test' ;;
esac

case $EXP5 in
Y|y)
# Bins Partition Sensitivity
nbinspar=0
for tmpbinspar in $binsparlst
do
  echo $tmpbinspar
  [[ -s idealcrtm.nml ]]&&rm idealcrtm.nml
  cat ${INNML} | sed s:_NAERS_:${naers}:g \
               | sed s:_USEREMI_:${useremi}:g \
               | sed s:_AERNAME_:${aername}:g \
               | sed s:_BINSLST_:${binlist}:g \
               | sed s:_BINSPAR_:${tmpbinspar}:g \
               | sed s:_GENMETHOD_:${genmethod}:g \
               | sed s:_TOTCONC_:${totconc}:g \
               | sed s:_THICKNS_:${thickns}:g \
               | sed s:_ATLEVEL_:${atlvl}:g \
               | sed s:_MODELV_:${modelv}:g \
               | sed s:_BMODERATIO_:${bmoderatio}:g \
               | sed s:_SFCMODRATIO_:${sfc_mode_ratio}:g \
               | sed s:_LANDCOV_:${landcover}:g \
               | sed s:_LANDTYPE_:${landtype}:g \
               | sed s:_LAI_:${lai}:g \
               | sed s:\'::g \
               > idealcrtm.nml
  if [ -s idealcrtm.nml ] ; then
     $EXEC
     rc=$?
     if [ $rc -eq 0 ]; then
        mkdir -p $OUTDIR/BinsPartition
        mv $KMODFILE $OUTDIR/BinsPartition/BinsPar${nbinspar}_kmod.nc
        mv $AODKFILE $OUTDIR/BinsPartition/BinsPar${nbinspar}_aodk.nc
     elif [ $rc -ne 0 ]; then
        exit 1
     fi
  fi
  let nbinspar=nbinspar+1
done
;;
*)
echo 'Skip Bins Partition test' ;;
esac

case $EXP6 in
Y|y)
# Surface type
for sfc in $sfctypelst
do
  case $sfc in
  'desert')
     explandcover='1.0' ; expsfctype='3'  ; explai='0.0' ;;
  'tundra')
     explandcover='1.0' ; expsfctype='10' ; explai='0.17';;
  'water')
     explandcover='0.0' ; expsfctype='10' ; explai='0.17';;
  *)
    echo 'not available surface type test' ;  exit ;;
  esac
  echo $sfc
  [[ -s idealcrtm.nml ]]&&rm idealcrtm.nml
  cat ${INNML} | sed s:_NAERS_:${naers}:g \
               | sed s:_USEREMI_:${useremi}:g \
               | sed s:_AERNAME_:${aername}:g \
               | sed s:_BINSLST_:${binlist}:g \
               | sed s:_BINSPAR_:${binspar}:g \
               | sed s:_GENMETHOD_:${genmethod}:g \
               | sed s:_TOTCONC_:${totconc}:g \
               | sed s:_THICKNS_:${thickns}:g \
               | sed s:_ATLEVEL_:${atlvl}:g \
               | sed s:_MODELV_:${modelv}:g \
               | sed s:_BMODERATIO_:${bmoderatio}:g \
               | sed s:_SFCMODRATIO_:${sfc_mode_ratio}:g \
               | sed s:_LANDCOV_:${explandcover}:g \
               | sed s:_LANDTYPE_:${expsfctype}:g \
               | sed s:_LAI_:${explai}:g \
               | sed s:\'::g \
               > idealcrtm.nml
  if [ -s idealcrtm.nml ] ; then
     $EXEC
     rc=$?
     if [ $rc -eq 0 ]; then
        mkdir -p $OUTDIR/SfcType
        mv $KMODFILE $OUTDIR/SfcType/${sfc}_kmod.nc
        mv $AODKFILE $OUTDIR/SfcType/${sfc}_aodk.nc
     elif [ $rc -ne 0 ]; then
        exit 1
     fi
  fi
done
;;
*)
echo 'Skip Surface Type test' ;;
esac

case $EXP7 in
Y|y)
  [[ -s idealcrtm.nml ]]&&rm idealcrtm.nml
  cat ${INNML} | sed s:_NAERS_:${naers}:g \
               | sed s:_USEREMI_:'.true.':g \
               | sed s:_AERNAME_:${aername}:g \
               | sed s:_BINSLST_:${binlist}:g \
               | sed s:_BINSPAR_:${binspar}:g \
               | sed s:_GENMETHOD_:${genmethod}:g \
               | sed s:_TOTCONC_:${totconc}:g \
               | sed s:_THICKNS_:${thickns}:g \
               | sed s:_ATLEVEL_:${atlevel}:g \
               | sed s:_MODELV_:${modelv}:g \
               | sed s:_BMODERATIO_:${bmoderatio}:g \
               | sed s:_SFCMODRATIO_:${sfc_mode_ratio}:g \
               | sed s:_LANDCOV_:${landcover}:g \
               | sed s:_LANDTYPE_:${landtype}:g \
               | sed s:_LAI_:${lai}:g \
               | sed s:_EMIMIN_:'0.8':g \
               | sed s:_EMIMAX_:'0.98':g \
               | sed s:_NEMISS_:'10':g \
               | sed s:\'::g \
               > idealcrtm.nml
  if [ -s idealcrtm.nml ] ; then
     $EXEC
     rc=$?
     if [ $rc -eq 0 ]; then
        mkdir -p $OUTDIR/UseEmis
        mv $KMODFILE $OUTDIR/UseEmis/UseEmis_kmod.nc
        mv $AODKFILE $OUTDIR/UseEmis/UseEmis_aodk.nc
     elif [ $rc -ne 0 ]; then
        exit 1
     fi
  fi
;;
*)
echo 'Skip Surface Emissivity test' ;;
esac
