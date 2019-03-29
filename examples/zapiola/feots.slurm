#!/bin/bash
#SBATCH --job-name=zapiola    # Job name
#SBATCH --ntasks=84                    # Run on a single CPU
#SBATCH --ntasks-per-node=84
#SBATCH --partition=highmem96
#SBATCH --time=8:00:00               # Time limit hrs:min:sec
#SBATCH --output=feots_%j.log   # Standard output and error log

module load gcc/7.4.0 openmpi/3.1.2/gcc/7.4.0 hdf5/1.10.3/parallel/openmpi/3.1.2/gcc/7.4.0 netcdf/4.6.1/parallel/openmpi/3.1.2/gcc/7.4.0

./GenMask
./RegionalExtraction
export OMP_NUM_THREADS=12
mpirun -np 7 -x OMP_NUM_THREADS ./FEOTSInitialize
mpirun -np 7 -x OMP_NUM_THREADS ./FEOTSDriver
