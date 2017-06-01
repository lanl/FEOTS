#!/bin/bash
#
#SBATCH --nodes=1 --nodefile=feotsnodes 
#SBATCH --job-name=feots_integrate
#SBATCH -joe
#

module purge
module load gcc openmpi

export OMP_NUM_THREADS=8

cd /home/${USER}/FEOTS/examples/gulfstream_dye/
date
mpirun -np 3 -x OMP_NUM_THREADS ./FEOTSDriver
date
