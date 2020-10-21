#!/bin/bash
echo "FEOTS_DBROOT=$FEOTS_DBROOT"
echo "REGIONAL_DB=$REGIONAL_DB"
echo "OUTDIR=$OUTDIR"

#bash ./bash/01-genmask.sh
#bash ./bash/02-region-extraction.sh

# Create this experiment OUTDIR
mkdir -p ${OUTDIR}
ln -s ${REGIONAL_DB}/* ${OUTDIR}/

bash ./bash/03-init.sh
bash ./bash/04-integrate.sh
bash ./bash/05-post-process.sh

echo "Mask Generation Job ID : $jid1"
echo "Region Extraction Job ID : $jid2"
echo "Initialization Job ID : $jid3"
echo "Integration Job ID : $jid4"
