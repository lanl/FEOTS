#!/bin/bash
echo "FEOTS_DBROOT=$FEOTS_DBROOT"
echo "REGIONAL_DB=$REGIONAL_DB"
echo "OUTDIR=$OUTDIR"

# Create this experiment OUTDIR
mkdir -p ${OUTDIR}
#bash ./bash/01-genmask.sh
#bash ./bash/02-genmaps.sh
#bash ./bash/03-init.sh
bash ./bash/04-equilibrate.sh
