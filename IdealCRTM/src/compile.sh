#!/bin/ksh

module purge
module load intel
module load netcdf
module list

BASE=/scratch1/BMC/gsd-fv3-dev/Shih-wei.Wei/IdealCRTM

export INCCRTMtest=$BASE/REL-2.3.0/crtm_v2.3.0_debug/crtm_v2.3.0/include
export LIBCRTMtest=$BASE/REL-2.3.0/crtm_v2.3.0_debug/crtm_v2.3.0/lib
#export INCCRTMtest=$BASE/REL-2.3.0/crtm_v2.3.0/crtm_v2.3.0/include
#export LIBCRTMtest=$BASE/REL-2.3.0/crtm_v2.3.0/crtm_v2.3.0/lib
export NETCDF_LIB=${NETCDF}/lib
export NETCDF_INC=${NETCDF}/include

make clean
make all

if [ -s crtm_ideal ]; then
   echo "compile succeed."
   mv crtm_ideal ../bin 
   make clean
else
   echo "compile error!!"
fi
