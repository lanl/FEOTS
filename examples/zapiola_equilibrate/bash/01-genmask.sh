#!/bin/bash

echo ""
echo "================================================"
echo "Creating mask for FEOTS simulation database"
echo ""
./genmask --dbroot ${FEOTS_DBROOT}  \
          --regional-db ${REGIONAL_DB} \
          --out ${OUTDIR} \
          --param-file ./runtime.params
