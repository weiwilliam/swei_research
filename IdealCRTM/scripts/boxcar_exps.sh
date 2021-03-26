#!/bin/ksh
#
HOMEDIR=/scratch1/BMC/gsd-fv3-dev/Shih-wei.Wei/IdealCRTM
USHDIR=${HOMEDIR}/ush
WRKDIR=${HOMEDIR}/wrk
OUTDIR=${HOMEDIR}/output
INNML=${USHDIR}/idealcrtm.nml.IN
EXEC=${HOMEDIR}/bin/crtm_ideal
KMODFILE='kmod_inspection.nc'
AODKFILE='aodk_inspection.nc'

# Exp on/off
EXP1='Y' # Binsize
EXP2='Y' # Concentration
EXP3='Y' # Thickness
EXP4='Y' # Altitude

# Module load
module load intel
module load netcdf

# Default setting for each factor
aername='DU002'
genmethod='2'

#genmethod=1
totconc='1.967e-03'
thickns='1'
atlevel='76'

daytime='.false.'

# genmethod=2
modelv='76'
extendlvs='4'
sfc_mode_ratio='0.5'

binlist="DU001 DU002 DU003 DU004 DU005"
conclist="9.835e-04 1.967e-03 3.934e-03 5.901e-03"
thicklist="1 3 5 7 9"
atlvllist="76 83 88 92"

rm -rf $WRKDIR/*
if [ ! -s $WRKDIR/coefficients ]; then
   mkdir -p $WRKDIR/coefficients
   cd $WRKDIR/coefficients
   sh $USHDIR/link_crtm_coeff.sh
fi
cd $WRKDIR

case $EXP1 in
Y|y)
# Bin Size Sensitivity
for aer in $binlist
do
  echo $aer
  [[ -s idealcrtm.nml ]]&&rm idealcrtm.nml
  cat ${INNML} | sed s:_AERNAME_:${aer}:g \
               | sed s:_GENMETHOD_:${genmethod}:g \
               | sed s:_TOTCONC_:${totconc}:g \
               | sed s:_THICKNS_:${thickns}:g \
               | sed s:_ATLEVEL_:${atlevel}:g \
               | sed s:_MODELV_:${modelv}:g \
               | sed s:_EXTENDLVS_:${extendlvs}:g \
               | sed s:_SFCMODRATIO_:${sfc_mode_ratio}:g \
               > idealcrtm.nml
              #| sed s:_DAYTIME_:${daytime}:g \
  if [ -s idealcrtm.nml ] ; then
     $EXEC
     rc=$?
     if [ $rc -eq 0 ]; then
        mkdir -p $OUTDIR/BinSize
        mv $KMODFILE $OUTDIR/BinSize/${aer}_kmod.nc
        mv $AODKFILE $OUTDIR/BinSize/${aer}_aodk.nc
     elif [ $rc -ne 0 ]; then
        exit 1
     fi
  fi
done
;;
*)
echo 'Skip bin size test' ;;
esac       

case $EXP2 in
Y|y)
# Total Concentration Sensitivity
for conc in $conclist
do
  echo $conc
  [[ -s idealcrtm.nml ]]&&rm idealcrtm.nml
  cat ${INNML} | sed s:_AERNAME_:${aername}:g \
               | sed s:_GENMETHOD_:${genmethod}:g \
               | sed s:_TOTCONC_:${conc}:g \
               | sed s:_THICKNS_:${thickns}:g \
               | sed s:_ATLEVEL_:${atlevel}:g \
               | sed s:_MODELV_:${modelv}:g \
               | sed s:_EXTENDLVS_:${extendlvs}:g \
               | sed s:_SFCMODRATIO_:${sfc_mode_ratio}:g \
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

case $EXP3 in
Y|y)
# Thickness Sensitivity
for thick in $thicklist
do
  echo $thick
  [[ -s idealcrtm.nml ]]&&rm idealcrtm.nml
  cat ${INNML} | sed s:_AERNAME_:${aername}:g \
               | sed s:_GENMETHOD_:${genmethod}:g \
               | sed s:_TOTCONC_:${totconc}:g \
               | sed s:_THICKNS_:${thick}:g \
               | sed s:_ATLEVEL_:${atlevel}:g \
               | sed s:_MODELV_:${modelv}:g \
               | sed s:_EXTENDLVS_:${extendlvs}:g \
               | sed s:_SFCMODRATIO_:${sfc_mode_ratio}:g \
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

case $EXP4 in
Y|y)
# Altitude Sensitivity
for atlvl in $atlvllist
do
  echo $atlvl
  [[ -s idealcrtm.nml ]]&&rm idealcrtm.nml
  cat ${INNML} | sed s:_AERNAME_:${aername}:g \
               | sed s:_GENMETHOD_:${genmethod}:g \
               | sed s:_TOTCONC_:${totconc}:g \
               | sed s:_THICKNS_:${thickns}:g \
               | sed s:_ATLEVEL_:${atlvl}:g \
               | sed s:_MODELV_:${modelv}:g \
               | sed s:_EXTENDLVS_:${extendlvs}:g \
               | sed s:_SFCMODRATIO_:${sfc_mode_ratio}:g \
               > idealcrtm.nml
  if [ -s idealcrtm.nml ] ; then
     $EXEC
     rc=$?
     if [ $rc -eq 0 ]; then
        mkdir -p $OUTDIR/Altitude
        mv $KMODFILE $OUTDIR/Altitude/L${thickns}_at_L${atlvl}_kmod.nc
        mv $AODKFILE $OUTDIR/Altitude/L${thickns}_at_L${atlvl}_aodk.nc
     elif [ $rc -ne 0 ]; then
        exit 1
     fi
  fi
done
;;
*)
echo 'Skip altitude test' ;;
esac
