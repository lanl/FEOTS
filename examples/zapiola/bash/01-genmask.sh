#!/bin/bash

#echo "================================================"
#echo "Creating mask for FEOTS regional database"
#echo ""
#feots genmask --dbroot ${FEOTS_DBROOT} \
#              --regional-db ${REGIONAL_DB} \
#              --param-file ./runtime.params
#
echo ""
echo "================================================"
echo "Creating mask for FEOTS simulation database"
echo ""
./genmask --dbroot ${FEOTS_DBROOT}  \
          --regional-db ${REGIONAL_DB} \
          --out ${OUTDIR} \
          --param-file ./runtime.params

