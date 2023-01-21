#!/bin/bash
#SBATCH --output=stdout
#SBATCH --error=stderr

#SBATCH --job-name=C14-b
#SBATCH --partition=std
#SBATCH --time=11:59:00

#SBATCH --nodes=6
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
initial_mode=51   # First phonon mode of this run
final_mode=102    # Last phonon mode of this run

###############################################################################################


export I_MPI_FABRICS=shm:ofa

rm -r Runs
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



# POSITIONS DISPLACED -

for j in $(seq $initial_mode $final_mode); do

   # Calculation of eigenvalues and wavefunctions
  mkdir mode-$j
  cp frozen-R/displaced-mode"$j"_scf.in  mode-$j
  cd mode-$j
  mpirun -np 96 $ESPRESSOLOCATION/pw.x < displaced-mode"$j"_scf.in > out_scf-mode"$j".out
  cp "$PREFIX".save/K00001/eigenval.xml ../outputfiles/outputsVeffcalc/displaced-/eigenval-mode"$j".xml
  mv out_scf-mode"$j".out ../outputfiles/outputsVeffcalc/displaced-

  # Printing of wavefunctions
  for k in $(seq $InitialBandFF $FinalBandFF); do
    # Modifying input files for every band
    cp ../inputfiles/wfcs.in ./wfcs-band$k.in
    i=' \ \ kband='$k','
    h=" \ \ filplot = '"$k".wfn',"
    sed -i "1i \ ${i}" wfcs-band$k.in
    sed -i "1i \ ${h}" wfcs-band$k.in
    sed -i '1i &InputPP' wfcs-band$k.in
    # Actual writing of the wavefunction
    mpirun -np 96 $ESPRESSOLOCATION/pp.x < wfcs-band$k.in > out_wfcs-band$k.out
    l=$((k-1))
    rm out_wfcs-band$l.out  # We remove these files because they are many and they demand a huge space
    rm wfcs-band$k.in
    mv $k.wfn mode"$j"-wfc"$k".wfn
    mv mode"$j"-wfc"$k".wfn ../Wfcs/Wfcs-displaced-/
  done

  cd ..
  rm -r mode-$j

done


