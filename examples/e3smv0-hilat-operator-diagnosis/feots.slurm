#!/bin/bash
#SBATCH --job-name=feots    # Job name
#SBATCH --ntasks=2                    # Run on a single CPU
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH --mem=150G
#SBATCH --time=1:00:00               # Time limit hrs:min:sec
#SBATCH --output=logs   # Standard output and error log
#SBATCH --account=feots
#SBATCH --partition=operator-diagnosis


module purge
module use /home/joe/modulefiles
module load gnu8/8.3.0 mpich/3.3.1 phdf5/1.10.5 netcdf/4.7.1 netcdf-fortran/4.5.2 feots/1.0.0


export OMP_NUM_THREADS=32
mkdir -p /mnt/feots/simulation/e3smv0-hilat-operator-diagnosis

#mpirun -n 2 -envall ./init --out /mnt/feots/simulation/e3smv0-hilat-operator-diagnosis

# Link the initial condition from NFS storage
#ln -s /mnt/feots/simulation/e3smv0-hilat-operator-diagnosis/Tracer.init.nc ./

mpirun -n 2 -envall feots integrate --dbroot /mnt/feots/db/E3SMV0-HILAT-5DAVG \
                                               --out /mnt/feots/simulation/e3smv0-hilat-operator-diagnosis
