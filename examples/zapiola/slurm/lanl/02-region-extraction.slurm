#!/bin/bash
#SBATCH --job-name=zapiola    # Job name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --time=1:00:00               # Time limit hrs:min:sec
#SBATCH --output=feots_logs   # Standard output and error log
#SBATCH -A climatehilat
#SBATCH --qos=standard
#SBATCH --array=365-365%1 # Modify the number of array indices to set how many operators to extract

module use /users/jschoonover/modulefiles
module load gcc/9.3.0 openmpi/3.1.6 hdf5-parallel/1.8.16 netcdf-h5parallel/4.4.0 feots/dev

feots region-extraction --dbroot ${FEOTS_DBROOT} \
                        --regional-db ${REGIONAL_DB} \
                        --oplevel ${SLURM_ARRAY_TASK_ID} \
                        --param-file ./runtime.params
