#!/bin/bash
#
# Uses the job-pipeline.sh to run time integration convergence tests for this example

set -x
dts=(1080)
#integrators=("euler")
integrators=("euler" "ab2" "ab3")
GCS_DEST="gs://feots-db/simulation/regional_constant_preservation"

export OUTDIR_ROOT=/home/joe/apps/feots_output/regional_constant_preservation/
export REGIONAL_DB=/home/joe/apps/feots_output/argentine_basin/
export FEOTS_DBROOT=/home/joe/apps/feots_output/E3SMV0-HILAT-5DAVG/

export FEOTS_FLAGS='--no-vertical-mixing'

# No vertical mixing tests
for dt in ${dts[@]}; do
  for integrator in ${integrators[@]}; do

    sed -i -e "s/timeStepScheme.*/timeStepScheme = '$integrator',/" ./runtime.params
    sed -i -e "s/dt.*/dt = $dt,/" ./runtime.params
    export OUTDIR="${OUTDIR_ROOT}/${integrator}-${dt}-nomixing"

    # Run the job-pipeline
    bash ./bash/job-pipeline.sh > ./feots.logs
    mv ./feots.logs ${OUTDIR}
    cp ./runtime.params ${OUTDIR}

    # Copy the output to a GCS bucket and remove the locally stored data (to save space)
    gsutil -m cp -r ${OUTDIR} ${GCS_DEST}
    rm -rf ${OUTDIR}

  done
done

## Vertical mixing tests
unset FEOTS_FLAGS
for dt in ${dts[@]}; do
  for integrator in ${integrators[@]}; do

    sed -i -e "s/timeStepScheme.*/timeStepScheme = '$integrator',/" ./runtime.params
    sed -i -e "s/dt.*/dt = $dt,/" ./runtime.params
    export OUTDIR="${OUTDIR_ROOT}/${integrator}-${dt}-mixing"
    bash ./bash/job-pipeline.sh > ./feots.logs
    mv ./feots.logs ${OUTDIR}

    gsutil cp -r ${OUTDIR} ${GCS_DEST}
    rm -rf ${OUTDIR}

  done
done
