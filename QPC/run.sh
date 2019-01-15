#!/bin/bash
set -ex
test -n "$1"
exe="$1" ; shift
sd=`pwd`

wd=/home2/ssd1/otkach/`date +%Y-%m-%d-2.5`
mkdir -p $wd && cd $wd

cp -p $sd/imp_*.dat .

cmd=(
  mpirun   
  -np 1 # number of processors (nodes?)
  -ppn 1 # processes per node
#  -t    # dry run
  -v    # verbose
  -maxtime 119  # minutes
  # -s broadwell
  $sd/$exe $*
)

#export OMP_NUM_THREADS=1
export KMP_AFFINITY=compact,granularity=fine

${cmd[*]}
