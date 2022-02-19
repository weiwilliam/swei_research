#!/bin/ksh
exefile=crtm_ideal_FDM
newexefile=crtm_ideal_FDM_v24
machine=$1
if [[ $machine == 'hera' ]]; then
   module purge
   module load intel
   module load netcdf
   module list
   BASE=/scratch1/BMC/gsd-fv3-dev/Shih-wei.Wei/IdealCRTM
   export INCCRTMtest=$BASE/REL-2.3.0/crtm_v2.3.0_debug/crtm_v2.3.0/include
   export LIBCRTMtest=$BASE/REL-2.3.0/crtm_v2.3.0_debug/crtm_v2.3.0/lib
   #export INCCRTMtest=$BASE/REL-2.3.0/crtm_v2.3.0/crtm_v2.3.0/include
   #export LIBCRTMtest=$BASE/REL-2.3.0/crtm_v2.3.0/crtm_v2.3.0/lib
elif [[ $machine == 's4' ]]; then
   module purge
   module load license_intel/S4
   module load intel
   module load hdf hdf5
   module load netcdf4
   module list
   BASE=/data/users/swei/Libs
   export INCCRTMtest=${BASE}/CRTM-2.4.0/include
   export LIBCRTMtest=${BASE}/CRTM-2.4.0/lib
   #export INCCRTMtest=${BASE}/REL-2.3.0/crtm_v2.3.0/crtm_v2.3.0/include
   #export LIBCRTMtest=${BASE}/REL-2.3.0/crtm_v2.3.0/crtm_v2.3.0/lib
   export NETCDF=$SSEC_NETCDF4_DIR
   expdir=/data/users/swei/Experiments/IdealCRTM
fi

export NETCDF_LIB=${NETCDF}/lib
export NETCDF_INC=${NETCDF}/include

make clean
make all

if [ -s $exefile ]; then
   echo "compile succeed."
   mv $exefile $expdir/bin/$newexefile 
   make clean
else
   echo "compile error!!"
fi
