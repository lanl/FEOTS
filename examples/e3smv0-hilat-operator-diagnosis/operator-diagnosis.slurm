#!/bin/bash
#SBATCH --job-name=op-diag    # Job name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --ntasks-per-node=1
#SBATCH --mem=210G
#SBATCH --time=1:00:00               # Time limit hrs:min:sec
#SBATCH --output=logs   # Standard output and error log
#SBATCH --account=feots
#SBATCH --partition=operator-diagnosis


module purge
module use /home/joe/modulefiles
module load gnu8/8.3.0 mpich/3.3.1 phdf5/1.10.5 netcdf/4.7.1 netcdf-fortran/4.5.2 feots/1.0.0


# Operator diagnosis
# Memory usage : ~201 GB, Single Core (Serial)
feots operator-diagnosis --dbroot /mnt/feots/db/E3SMV0-HILAT-5DAVG \
                         --irf /mnt/feots/db/E3SMV0-HILAT-5DAVG/irf/response/IRF.0064-01-05.nc \
                         --oplevel 1 
