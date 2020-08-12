#!/bin/bash
#SBATCH --job-name=op-diag    # Job name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=215G
#SBATCH --time=1:00:00               # Time limit hrs:min:sec
#SBATCH --output=logs-%j.txt   # Standard output and error log
#SBATCH --account=feots
#SBATCH --partition=operator-diagnosis
###SBATCH --array=199-219%1
#SBATCH --array=3-4%1


module purge
module use /home/joe/modulefiles
module load gnu8/8.3.0 mpich/3.3.1 phdf5/1.10.5 netcdf/4.7.1 netcdf-fortran/4.5.2 feots/1.0.0

GCSROOT=gs://feots-db/db/E3SMV0-HILAT-5DAVG/
DBROOT=/mnt/feots/db/E3SMV0-HILAT-5DAVG/
IRF=$(sed -n ${SLURM_ARRAY_TASK_ID}p IRFs.txt)
echo $IRF
# Operator diagnosis
# Memory usage : ~201 GB, Single Core (Serial)
feots operator-diagnosis --dbroot $DBROOT \
                         --irf $IRF \
                         --oplevel $SLURM_ARRAY_TASK_ID 

# Copy the resultant operators to GCS Bucket
OPID=$(printf "%05d" $SLURM_ARRAY_TASK_ID)
gsutil cp ${DBROOT}/ops/diffusion.${OPID}.h5 ${GCSROOT}/ops/ 
gsutil cp ${DBROOT}/ops/transport.${OPID}.h5 ${GCSROOT}/ops/ 
