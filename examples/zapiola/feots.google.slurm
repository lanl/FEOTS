#!/bin/bash
#SBATCH --job-name=regional-extraction    # Job name
#SBATCH --ntasks=7                    # Run on a single CPU
#SBATCH --ntasks-per-node=7
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --time=1:00:00               # Time limit hrs:min:sec
#SBATCH --output=logs   # Standard output and error log
#SBATCH --account=feots
#SBATCH --partition=operator-diagnosis

module purge
module use /home/joe/modulefiles
module load gnu8/8.3.0 mpich/3.3.1 phdf5/1.10.5 netcdf/4.7.1 netcdf-fortran/4.5.2 feots/1.0.0

#mpirun -n 7 ./init --dbroot /mnt/feots/db/E3SMV0-HILAT-5DAVG \
#                   --param-file ./runtime.google.params \
#                   --out /mnt/feots/solution/zapiola

ln -s /mnt/feots/solution/zapiola/Tracer.init.nc ./

mpirun -n 7 feots integrate --dbroot /mnt/feots/db/E3SMV0-HILAT-5DAVG \
                            --param-file ./runtime.google.params \
                            --out /mnt/feots/solution/zapiola
