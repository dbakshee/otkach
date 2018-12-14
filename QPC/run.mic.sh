#!/bin/bash
set -ex
test -n "$1"
exe="$1" ; shift
dirname=`date +%Y-%m-%d`
mkdir -p $dirname && cd $dirname

cp ../imp_*.dat .
echo node1:0  >hosts
echo node2:1 >>hosts
cmd=(
  mpirun 
  -np 1 # number of processors (nodes?)
#  -ppn 1 # processes per node
#  -t    # dry run
  -v    # verbose
  -termtime 1   # minutes
  -maxtime 10   # minutes
  -mic
  -machinefile hosts
  -mic-postfix .mic
  -s omnipath
  ../$exe $*
)
#export OMP_NUM_THREADS=16
export KMP_AFFINITY=compact,granularity=fine

"${cmd[@]}"

