#!/bin/bash
# LAUNCH THIS FILE WITH: mpirun -np 64 -ppn 64 -maxtime 1 ./run.sh ...
# usage: dis-bar.py [-h] [--vh VH] [--npi NPI]
#                   --temperature TEMPERATURE --E0 E0 --nkt NKT --ne NE --V0 BEG,END,N
#                   [--vr VR] [--Bt BT]
. ~/.bashrc
conda3
cd /home2/ifpso1/otkach/MagBreak
set -ex
export OMP_NUM_THREADS=1
chmod +x MagBreak.py
env | grep RANK

mkdir -p $SUPPZ_JOB_NAME
if [ "$PMI_RANK" == 0 ]; then
    cp $0 ./MagBreak.py $SUPPZ_JOB_NAME/
fi
cd $SUPPZ_JOB_NAME

for E in $* ; do
   ../MagBreak.py -E $E -I 5,35,1001 -vr 1 -w0 0.5
done
