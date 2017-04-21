#!/bin/bash
#
#SBATCH --nodes=1 --nodefile=feotsnodes 
#SBATCH --job-name=feots_integrate
#SBATCH -joe
#

module purge
module load gcc

cd /home/schoonover/feots_lite/examples/agulhas/
date
../../bin/FEOTSDriver
date
