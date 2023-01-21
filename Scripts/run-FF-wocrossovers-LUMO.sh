#!/bin/bash
#SBATCH --output=stdout
#SBATCH --error=stderr

#SBATCH --job-name=C14LUMO
#SBATCH --partition=std
#SBATCH --time=05:59:00

#SBATCH --nodes=1
#SBATCH --tasks-per-node=1

##you need to execute the .sh in the directory /work/fcnv733 ($WORK)

#set -e ###This stops when there is an error; also if you try to create an already existing directory

module switch env env/intel-15.0.3_impi-5.0.3

ifort -o frozen-phononLUMO.x frozen-phonon.f90

./frozen-phononLUMO.x  > out-FF-LUMO-wocrossovers.out

