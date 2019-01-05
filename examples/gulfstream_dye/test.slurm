#!/bin/bash
#SBATCH --job-name=feots    # Job name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --ntasks-per-node=1
#SBATCH --partition=highmem96
#SBATCH --time=2:00:00               # Time limit hrs:min:sec
#SBATCH --output=feots_%j.log   # Standard output and error log

module load gcc/7.4.0 openmpi/3.1.2/gcc/7.4.0 hdf5/1.10.3/parallel/openmpi/3.1.2/gcc/7.4.0 netcdf/4.6.1/parallel/openmpi/3.1.2/gcc/7.4.0

#./GenerateMeshOnlyFile
#./GenMask
#./OperatorDiagnosis
./RegionalExtraction



