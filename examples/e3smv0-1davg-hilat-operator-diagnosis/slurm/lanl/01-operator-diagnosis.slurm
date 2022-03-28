#!/bin/bash
#SBATCH --job-name=op-diag    # Job name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=56G
#SBATCH --time=15:00               # Time limit hrs:min:sec
#SBATCH --output=feots_logs-%j   # Standard output and error log
#SBATCH --account=climatehilat
#SBATCH --qos=standard
#SBATCH --array=2-105%2

module purge
module use /users/jschoonover/modulefiles
module load gcc/8.3.0 mvapich2/2.3 hdf5-parallel/1.8.16 netcdf-h5parallel/4.4.0 feots/dev

echo "Set FEOTS Database root to ${FEOTS_DBROOT}"
IRF=$(sed -n ${SLURM_ARRAY_TASK_ID}p IRFs.txt)
echo $IRF
feots operator-diagnosis --dbroot $FEOTS_DBROOT \
                         --irf $IRF \
                         --oplevel $SLURM_ARRAY_TASK_ID 

