#!/bin/bash
#SBATCH --job-name=zapiola    # Job name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --array=133-365
#SBATCH --time=1:00:00               # Time limit hrs:min:sec
#SBATCH --output=feots_%j.log   # Standard output and error log
#SBATCH -A climatehilat
#SBATCH --qos=standard

module use /users/jschoonover/modulefiles
module load gcc/8.3.0 mvapich2/2.3 hdf5-parallel/1.8.16 netcdf-h5parallel/4.4.0 feots/1.0.0


./process_files.sh runtime.params IRFList_5dayAvg.txt ${SLURM_ARRAY_TASK_ID}


#./genmask
#feots region-extraction
#export OMP_NUM_THREADS=12
#mpirun -np 7 -x OMP_NUM_THREADS ./init
#mpirun -np 7 -x OMP_NUM_THREADS feots integrate
#mpirun -np 7 -x OMP_NUM_THREADS feots equilibrate
