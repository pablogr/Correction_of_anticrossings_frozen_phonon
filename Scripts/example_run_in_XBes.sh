#!/bin/bash

#$ -e stderr
#$ -o stdout

#$ -cwd
#$ -pe parallel 32 # number of processors
#$ -q main.q      # name of queue
#$ -N C14cut30 # name of the process_programm which could be seen in the cluster
#$ -M pablo.risueno@chemie.uni-hamburg.de #your email for system reports

export I_MPI_FABRICS=shm:ofa

rm -r Runs
mkdir Wfcs
mkdir inputfiles
mkdir inputfiles/bandwritinginputfiles
mkdir outputfiles
mkdir Runs
mkdir Runs/dyneq


# STEP1: GEOMETRY OPTIMIZATION
mpirun -np 32 /home/risueno/CalcsPengProject/QE_files/pw.x < inputfiles/geomopt.in > outputfiles/out_geomopt.out  

echo ' STEP 1 (GEOMETRY OPTIMIZATION) FINISHED'
