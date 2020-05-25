#!/bin/bash
#SBATCH --job-name=zapiola    # Job name
#SBATCH --ntasks=84                    # Run on a single CPU
#SBATCH --ntasks-per-node=84
#SBATCH --partition=highmem96
#SBATCH --time=8:00:00               # Time limit hrs:min:sec
#SBATCH --output=feots_%j.log   # Standard output and error log

module use /users/jschoonover/feots
module load gcc/8.3.0 mvapich2/2.3 hdf5-parallel/1.8.16 netcdf-h5parallel/4.4.0 feots/1.0.0

./genmask
feots region-extraction
export OMP_NUM_THREADS=12
mpirun -np 7 -x OMP_NUM_THREADS ./init
mpirun -np 7 -x OMP_NUM_THREADS feots integrate
#mpirun -np 7 -x OMP_NUM_THREADS feots equilibrate
