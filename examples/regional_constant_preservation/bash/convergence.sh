#!/bin/bash
#
# Uses the job-pipeline.sh to run time integration convergence tests for this example

dts=(1080 540 270 135)
integrators=("euler" "ab2" "ab3")

export REGIONAL_DB=/home/joe/apps/feots_output/regional_constant_preservation/
export FEOTS_FLAGS='--no-vertical-mixing'

# No vertical mixing tests
for dt in ${dts[@]}; do
  for integrator in ${integrators[@]}; do

    sed -i -e "s/timeStepScheme.*/timeStepScheme = '$integrator',/" ./runtime.params
    sed -i -e "s/dt.*/dt = $dt,/" ./runtime.params
    export OUTDIR="${REGIONAL_DB}/${integrator}-${dt}-nomixing"
    bash ./bash/job-pipeline.sh > ./feots.logs
    mv ./feots.logs ${OUTDIR}

  done
done

## Vertical mixing tests
#unset FEOTS_FLAGS
#for dt in ${dts[@]}; do
#  for integrator in ${integrators[@]}; do
#
#    sed -i -e "s/timeStepScheme.*/timeStepScheme = '$integrator',/" ./runtime.params
#    sed -i -e "s/dt.*/dt = $dt,/" ./runtime.params
#    export OUTDIR="${REGIONAL_DB}/${integrator}-${dt}-mixing"
#    bash ./bash/job-pipeline.sh > ./feots.logs
#    mv ./feots.logs ${OUTDIR}
#
#  done
#done
