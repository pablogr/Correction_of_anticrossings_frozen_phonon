#!/bin/bash
#SBATCH --output=stdout
#SBATCH --error=stderr

#SBATCH --job-name=goC14c30
#SBATCH --partition=std
#SBATCH --time=11:00:00

#SBATCH --nodes=6
#SBATCH --tasks-per-node=16

##you need to execute the .sh in the directory /work/fcnv733 ($WORK)

#set -e ###This stops when there is an error; also if you try to create an already existing directory
set -x

. /sw/batch/init.sh

result=$PWD
workdir=$RRZ_GLOBAL_TMPDIR

export ESPRESSOLOCATION=/home/fcnv733/PROGRAMS/espressoHahn/bin/

export LD_LIBRARY_PATH=/home/fcnv733/PROGRAMS/libmkl_PGR/


module switch env env/intel-15.0.3_impi-5.0.3


export I_MPI_FABRICS=shm:ofa

rm -r Runs
mkdir Wfcs
mkdir inputfiles
mkdir inputfiles/bandwritinginputfiles
mkdir outputfiles
mkdir Runs
mkdir Runs/Veff
mkdir Runs/dyneq


# STEP1: GEOMETRY OPTIMIZATION

mpirun -np 96 $ESPRESSOLOCATION/pw.x < inputfiles/geomopt.in > outputfiles/out_geomopt.out  

echo ' STEP 1 (GEOMETRY OPTIMIZATION) FINISHED'

