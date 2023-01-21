#!/bin/bash
#SBATCH --output=stdout
#SBATCH --error=stderr

#SBATCH --job-name=C14wfc0
#SBATCH --partition=std
#SBATCH --time=6:00:00

#SBATCH --nodes=2
#SBATCH --tasks-per-node=16

##you need to execute the .sh in the directory /work/fcnv733 ($WORK)

#set -e ###This stops when there is an error; also if you try to create an already existing directory
set -x

. /sw/batch/init.sh

result=$PWD
workdir=$RRZ_GLOBAL_TMPDIR


### Variables and paths (change their value to the value appropriate for your system) #######

export LD_LIBRARY_PATH=/home/fcnv733/PROGRAMS/libmkl_PGR/
module switch env env/intel-15.0.3_impi-5.0.3
export ESPRESSOLOCATION=/home/fcnv733/PROGRAMS/espressoHahn/bin/
export PREFIX=diam-

InitialBandFF=36 # Initial electronic orbital with overlap with HOMO and/or LUMO is to be checked
FinalBandFF=39   # Final electronic orbital whose overlap with HOMO and/or LUMO is to be checked

###############################################################################################


export I_MPI_FABRICS=shm:ofa

rm -r Runs
mkdir Veffs
mkdir Wfcs
mkdir Wfcs/Wfcs-undisplaced
mkdir Wfcs/Wfcs-displaced+
mkdir Wfcs/Wfcs-displaced-
mkdir inputfiles
mkdir inputfiles/bandwritinginputfiles
mkdir outputfiles
mkdir outputfiles/outputsVeffcalc
mkdir outputfiles/outputsVeffcalc/displaced+
mkdir outputfiles/outputsVeffcalc/displaced-
mkdir Runs
mkdir Runs/dyneq


# UNDISPLACED POSITIONS
mpirun -np 32 $ESPRESSOLOCATION/pw.x < inputfiles/scf.in > outputfiles/out_scf0.out
cp "$PREFIX".save/K00001/eigenval.xml .
cp outputfiles/out_scf0.out .

# Writing input files
for j in $(seq $InitialBandFF $FinalBandFF); do
  # Modifying input files for every band
  cp inputfiles/wfcs.in wfcs-band$j.in
  i=' \ \ kband='$j','
  k=" \ \ filplot = '"$j".wfn'," 
  sed -i "1i \ ${i}" wfcs-band$j.in
  sed -i "1i \ ${k}" wfcs-band$j.in
  sed -i '1i &InputPP' wfcs-band$j.in
  # Actual writing of the wavefunction
  mpirun -np 32 $ESPRESSOLOCATION/pp.x < wfcs-band$j.in > outputfiles/out_wfcs-band$j.out
  l=$((j-1))
  rm outputfiles/out_wfcs-band$l.out  # We remove these files because they are many and they demand a huge space
  mv wfcs-band$j.in inputfiles/bandwritinginputfiles
  mv $j.wfn wfc$j.wfn
  mv wfc$j.wfn Wfcs/Wfcs-undisplaced
done  


