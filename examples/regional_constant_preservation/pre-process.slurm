#!/bin/bash
#SBATCH --job-name=zapiola    # Job name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --cpus-per-task=1
#SBATCH --mem=30G
#SBATCH --array=1-5
#SBATCH --time=1:00:00               # Time limit hrs:min:sec
#SBATCH --output=feots_%j.log   # Standard output and error log
#SBATCH -A climatehilat
#SBATCH --qos=standard

module use /users/jschoonover/modulefiles
module load gcc/8.3.0 mvapich2/2.3 hdf5-parallel/1.8.16 netcdf-h5parallel/4.4.0 feots/1.0.0


./genmask
./process_files.sh runtime.params IRFList_5dayAvg.txt ${SLURM_ARRAY_TASK_ID}

