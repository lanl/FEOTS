#!/bin/bash
#SBATCH --job-name=zapiola    # Job name
#SBATCH --ntasks=6                    # Run on a single CPU
#SBATCH --cpus-per-task=1
#SBATCH --mem=30G
#SBATCH --time=1:00:00               # Time limit hrs:min:sec
#SBATCH --output=feots_logs   # Standard output and error log
#SBATCH -A climatehilat
#SBATCH --qos=standard

module use /users/jschoonover/modulefiles
module load gcc/9.3.0 openmpi/3.1.6 hdf5-parallel/1.8.16 netcdf-h5parallel/4.4.0 feots/dev

echo $SLURM_NTASKS
mpirun -np ${SLURM_NTASKS} ./init --out ${OUTDIR} \
                    --regional-db ${REGIONAL_DB} \
                    --dbroot ${FEOTS_DBROOT} \
                    --param-file ./runtime.params
