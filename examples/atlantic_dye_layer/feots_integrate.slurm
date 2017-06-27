#!/bin/bash
#
#SBATCH --nodes=1 --nodefile=feotsnodes 
#SBATCH --job-name=feots_integrate
#SBATCH -joe
#SBATCH --time=40:00:00
#

module purge
module load gcc openmpi

export OMP_NUM_THREADS=3

cd /home/${USER}/FEOTS/examples/atlantic_dye/
date
mpirun -np 4 -x OMP_NUM_THREADS ./FEOTSDriver	
				# np = number of processors
				# = # of tracers + 1
date
