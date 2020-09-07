#!/bin/bash
#SBATCH --job-name=op-diag    # Job name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=85G
#SBATCH --time=15:00               # Time limit hrs:min:sec
#SBATCH --output=feots_logs   # Standard output and error log
#SBATCH --account=climatehilat
#SBATCH --qos=standard
#SBATCH --array=2-5%1


module purge
module use /users/jschoonover/modulefiles
module load gcc/8.3.0 mvapich2/2.3 hdf5-parallel/1.8.16 netcdf-h5parallel/4.4.0 feots/dev

#GCSROOT=gs://feots-db/db/E3SMV0-HILAT-5DAVG/
#DBROOT=/mnt/feots/db/E3SMV0-HILAT-5DAVG/
DBROOT=/lustre/scratch3/turquoise/jschoonover/feots/E3SMV0-HILAT-5DAVG
OUTDIR=/lustre/scratch3/turquoise/${USER}/global

IRF=$(sed -n ${SLURM_ARRAY_TASK_ID}p IRFs.txt)
echo $IRF
# Operator diagnosis
# Memory usage : ~201 GB, Single Core (Serial)
feots operator-diagnosis --dbroot $DBROOT \
                         --irf $IRF \
                         --oplevel $SLURM_ARRAY_TASK_ID 

