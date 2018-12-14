#!/bin/bash
set -ex
test -n "$1"
exe="$1" ; shift
dirname=`date +%Y-%m-%da`
mkdir -p $dirname && cd $dirname
ln -f -s ../imp_*.dat .
cmd=(
  mpirun 
  -np 1 # number of processors (nodes?)
  -ppn 1 # processes per node
#  -t    # dry run
  -v    # verbose
  -maxtime 300   # minutes
  ../$exe $*
)
#export OMP_NUM_THREADS=16
export KMP_AFFINITY=compact,granularity=fine

${cmd[*]}
