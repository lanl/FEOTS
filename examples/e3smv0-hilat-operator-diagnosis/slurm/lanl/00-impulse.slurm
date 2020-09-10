#!/bin/bash
#SBATCH --job-name=op-diag    # Job name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=56G
#SBATCH --time=5:00               # Time limit hrs:min:sec
#SBATCH --output=feots_logs   # Standard output and error log
#SBATCH --account=climatehilat
#SBATCH --qos=standard

module purge
module use /users/jschoonover/modulefiles
module load gcc/8.3.0 mvapich2/2.3 hdf5-parallel/1.8.16 netcdf-h5parallel/4.4.0 feots/dev

#DBROOT=/lustre/scratch3/turquoise/jschoonover/feots/E3SMV0-HILAT-5DAVG

feots impulse --dbroot $FEOTS_DBROOT \
              --irf $POP_MESH_FILE

