#!/bin/ksh

CRTMfix=/scratch1/BMC/gsd-fv3-dev/Shih-wei.Wei/IdealCRTM/REL-2.3.0/fix

sensorlist="iasi_metop-a v.modis_aqua"
BoL="Big" # Big, Little
AbsAlgorithm='ODPS'

SLN="ln -fs"

$SLN ${CRTMfix}/AerosolCoeff/${BoL}_Endian/* .
$SLN ${CRTMfix}/CloudCoeff/${BoL}_Endian/* .
$SLN ${CRTMfix}/EmisCoeff/IR_Ice/SEcategory/${BoL}_Endian/* .
$SLN ${CRTMfix}/EmisCoeff/IR_Land/SEcategory/${BoL}_Endian/* .
$SLN ${CRTMfix}/EmisCoeff/IR_Snow/SEcategory/${BoL}_Endian/* .
$SLN ${CRTMfix}/EmisCoeff/IR_Water/${BoL}_Endian/* .
$SLN ${CRTMfix}/EmisCoeff/MW_Water/${BoL}_Endian/* .
$SLN ${CRTMfix}/EmisCoeff/VIS*/SEcategory/${BoL}_Endian/* .
for sensor in $sensorlist
do
  case $sensor in
  "iasi_metop-a")
     prefix="iasi_metop-a" ;;
  *)
     prefix=$sensor
  esac
  $SLN ${CRTMfix}/SpcCoeff/${BoL}_Endian/${prefix}.SpcCoeff.bin .
  if [ -s ${CRTMfix}/TauCoeff/${AbsAlgorithm}/${BoL}_Endian/${prefix}.TauCoeff.bin ]; then
     $SLN ${CRTMfix}/TauCoeff/${AbsAlgorithm}/${BoL}_Endian/${prefix}.TauCoeff.bin .
  else
     $SLN ${CRTMfix}/TauCoeff/ODAS/${BoL}_Endian/${prefix}.TauCoeff.bin .
  fi
done
