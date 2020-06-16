#!/bin/bash
#SBATCH --job-name=constant_preservation    # Job name
#SBATCH --ntasks=2                    # Run on a single CPU
#SBATCH --cpus-per-task=12
#SBATCH --mem=120G
#SBATCH --time=1:00:00               # Time limit hrs:min:sec
#SBATCH --output=feots_%j.log   # Standard output and error log
#SBATCH -A climatehilat
#SBATCH --qos=standard

module purge
module use /users/jschoonover/modulefiles
#module load gcc/8.3.0 mvapich2/2.3 hdf5-parallel/1.8.16 netcdf-h5parallel/4.4.0 feots/1.0.0
module load gcc/8.3.0 mvapich2/2.3 hdf5-parallel/1.8.16 netcdf-h5parallel/4.4.0 feots/dev

#module load gcc/8.3.0 hdf5-serial/1.8.16 netcdf-serial/4.4.0 feots/gcc-8.3.0/1.0.0

export OMP_NUM_THREADS=12
#srun -n2 -x OMP_NUM_THREADS ./init
srun -n2 -x OMP_NUM_THREADS feots integrate
#srun -n2 ./init
#srun -n2 feots integrate

#./init
