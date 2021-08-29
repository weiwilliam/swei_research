#!/bin/ksh

export machine='s4'

if [ $machine == 'hera' ] ; then
   gsisub=sub_hera
   account=star
   gsirunscript=$1
   procs=240/1
   wtime='01:00:00'
   outfile=/scratch1/BMC/gsd-fv3-dev/Shih-wei.Wei/runlogs/log.$gsirunscript
elif [ $machine == 's4' ] ; then
   gsisub=/data/users/swei/Git/GSI/ush/sub_s4
   account=star
   gsirunscript=$1
   # nprocs per node / nodes
   procs='12/12'
   #procs='12/16'
   wtime='00:30:00'
   outfile=/data/users/swei/Experiments/runlogs/log.$gsirunscript
fi

$gsisub -j $gsirunscript -o $outfile -a $account -p $procs -t $wtime $gsirunscript


