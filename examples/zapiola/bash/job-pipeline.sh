#!/bin/bash
echo "FEOTS_DBROOT=$FEOTS_DBROOT"
echo "REGIONAL_DB=$REGIONAL_DB"
echo "OUTDIR=$OUTDIR"

#bash ./bash/01-genmask.sh
#bash ./bash/02-region-extraction.sh

# Create this experiment OUTDIR
mkdir -p ${OUTDIR}

set -x
bash ./bash/01-genmask.sh
bash ./bash/02-genmaps.sh
bash ./bash/03-init.sh
bash ./bash/04-integrate.sh
#bash ./bash/05-post-process.sh
