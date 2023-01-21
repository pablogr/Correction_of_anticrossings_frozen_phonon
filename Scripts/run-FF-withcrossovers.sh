#!/bin/bash
#SBATCH --output=stdout
#SBATCH --error=stderr

#SBATCH --job-name=C14FF
#SBATCH --partition=std
#SBATCH --time=11:59:00

#SBATCH --nodes=1
#SBATCH --tasks-per-node=16

##you need to execute the .sh in the directory /work/fcnv733 ($WORK)


module switch env env/intel-15.0.3_impi-5.0.3

ifort  -o frozen-phonon.x frozen-phonon.f90
#ifort -g -traceback -o frozen-phonon.x frozen-phonon.f90

./frozen-phonon.x  > out-FF-withcrossovers.out
