#!/bin/bash
#
#SBATCH --nodes=1 --nodefile=feotsnodes 
#SBATCH --job-name=feots_diagnosis
#SBATCH -joe
#

module purge
module load gcc

cd /home/${USER}/FEOTS/examples/global_dye/
date
./OperatorDiagnosis
date
