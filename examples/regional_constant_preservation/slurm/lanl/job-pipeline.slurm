#!/bin/bash
echo "FEOTS_DBROOT=$FEOTS_DBROOT"
echo "OUTDIR=$OUTDIR"

#jid1=$(sbatch --get-user-env slurm/lanl/01-genmask.slurm | awk '{print $4}')
jid2=$(sbatch --get-user-env --dependency=afterok:$jid1 slurm/lanl/02-region-extraction.slurm | awk '{print $4}')
jid3=$(sbatch --get-user-env --dependency=afterok:$jid2 slurm/lanl/03-init.slurm | awk '{print $4}')
jid4=$(sbatch --get-user-env --dependency=afterok:$jid3 slurm/lanl/04-integrate.slurm | awk '{print $4}')

echo "Mask Generation Job ID : $jid1"
echo "Region Extraction Job ID : $jid2"
echo "Initialization Job ID : $jid3"
echo "Integration Job ID : $jid4"
