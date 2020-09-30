#!/bin/bash
# request resources:
#PBS -N vcisek-benchmark
#PBS -l nodes=1:ppn=8
#PBS -l walltime=00:05:00


if [ $PBS_O_WORKDIR ]; then
    cd $PBS_O_WORKDIR
fi

python simulate.py

