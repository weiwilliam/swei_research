#!/bin/ksh

JOBSQUEUE="`which squeue` -u ${USER}"
SQFORMAT="%.10i %.9P %.25j %.8u %.8T %.10M %.10L %.3D %R"
ndatepy=$HOME/bin/ndate.py

export machine='s4'

if [ $machine == 'hera' ] ; then
   gsisub=/scratch2/BMC/gsd-fv3-dev/Shih-wei.Wei/GSI/ush/sub_hera
   account=gsd-fv3-dev
   queue=debug
   gsirunscript=$1
   procs=20/8
   wtime='00:30:00'
   outfile="/scratch1/BMC/gsd-fv3-dev/Shih-wei.Wei/wrklog/${gsirunscript}.runlog.%j"
elif [ $machine == 's4' ] ; then
   gsisub=/data/users/swei/Git/GSI/ush/sub_s4
   account=star
   queue=debug
   #queue=batch
   gsirunscript=$1
   # nprocs per node / nodes
   procs='20/8'
   #procs='12/16'
   wtime='00:15:00'
   outfile=/data/users/swei/Experiments/runlogs/log.$gsirunscript
fi

SDATE=2020060112
EDATE=2020060112
CDATE=$SDATE

until [ $CDATE -gt $EDATE ]; do
   export CDATE=$CDATE
   $gsisub -j $gsirunscript -o $outfile -a $account -p $procs -t $wtime -q $queue $gsirunscript
   CDATE=`python $ndatepy 6 $CDATE`
   sqrc=0
   until [ $sqrc -ne 0 ]; do
      $JOBSQUEUE -o "${SQFORMAT}" | grep "$gsirunscript"
      sqrc=$?
      sleep 10
   done
done

