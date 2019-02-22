#!/bin/bash
#SBATCH --job-name=feots    # Job name
#SBATCH --ntasks=24                    # Run on a single CPU
#SBATCH --ntasks-per-node=24
#SBATCH --partition=highmem96
###SBATCH --time=8:00:00               # Time limit hrs:min:sec
#SBATCH --output=feots_%j.log   # Standard output and error log

module load gcc/7.4.0 openmpi/3.1.2/gcc/7.4.0 hdf5/1.10.3/parallel/openmpi/3.1.2/gcc/7.4.0 netcdf/4.6.1/parallel/openmpi/3.1.2/gcc/7.4.0

#./GenerateMeshOnlyFile

#./GenMask
./OperatorDiagnosis
./RegionalExtraction
#

export OMP_NUM_THREADS=20
mpirun -np 3 -x OMP_NUM_THREADS ./FEOTSInitialize
mpirun -np 3 -x OMP_NUM_THREADS ./FEOTSDriver


GIT_HASH=$(git log --pretty=format:'%h' -n 1)

mkdir -p hash_${GIT_HASH}/forward_input-decks
mkdir -p hash_${GIT_HASH}/forward_netcdf
mkdir -p hash_${GIT_HASH}/plots

cp *_mask.nc hash_${GIT_HASH}/forward_input-decks
cp *_mesh.nc hash_${GIT_HASH}/forward_input-decks
cp runtime.params hash_${GIT_HASH}/forward_input-decks

mv Tracer.*.nc hash_${GIT_HASH}/forward_netcdf

