#!/bin/bash
#SBATCH --output=stdout
#SBATCH --error=stderr

#SBATCH --job-name=dynC14Cc30
#SBATCH --partition=std
#SBATCH --time=11:59:00

#SBATCH --nodes=8
#SBATCH --tasks-per-node=16

##you need to execute the .sh in the directory /work/fcnv733 ($WORK)

#set -e ###This stops when there is an error; also if you try to create an already existing directory
set -x

. /sw/batch/init.sh

result=$PWD
workdir=$RRZ_GLOBAL_TMPDIR

export LD_LIBRARY_PATH=/home/fcnv733/PROGRAMS/libmkl_PGR/
export ESPRESSOLOCATION=/home/fcnv733/PROGRAMS/espressoHahn/bin/

module switch env env/intel-15.0.3_impi-5.0.3


rm -r Runs
mkdir Wfcs
mkdir inputfiles
mkdir inputfiles/bandwritinginputfiles
mkdir outputfiles
mkdir Runs
mkdir Runs/dyneq


mpirun -np 128 $ESPRESSOLOCATION/pw.x < inputfiles/scf.in > outputfiles/out_scf.out
mpirun -np 128 $ESPRESSOLOCATION/ph.x < inputfiles/dyneq.in > outputfiles/out_dyneq.out
mpirun -np 128 $ESPRESSOLOCATION/dynmat.x < inputfiles/dynmat.in > outputfiles/out_dynmat.out

