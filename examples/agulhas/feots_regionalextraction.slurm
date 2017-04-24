#!/bin/bash
#
#SBATCH --nodes=1 --nodefile=feotsnodes 
#SBATCH --job-name=feots_extraction
#SBATCH -joe
#

module purge
module load gcc

cd /home/${USER}/FEOTS/examples/agulhas/
date
./RegionalExtraction
date
