#!/bin/bash
#SBATCH --job-name=zapiola    # Job name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --time=1:00:00               # Time limit hrs:min:sec
#SBATCH --output=feots_logs   # Standard output and error log
#SBATCH -A climatehilat
#SBATCH --qos=standard
#SBATCH --array=1-100%20 # Modify the number of array indices to set how many operators to extract

module use /users/jschoonover/modulefiles
module load gcc/8.3.0 mvapich2/2.3 hdf5-parallel/1.8.16 netcdf-h5parallel/4.4.0 feots/opt

feots region-extraction --dbroot ${FEOTS_DBROOT} \
                        --out ${OUTDIR} \
                        --oplevel ${SLURM_ARRAY_TASK_ID} \
                        --param-file ./params/lanl/runtime.params
