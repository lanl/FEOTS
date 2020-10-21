# FEOTS - Zapiola Demo

## Description



## Directory Structure
This example contains directories for 
1. `params/` - Contains namelist parameter files for driving FEOTS simulations
2. `slurm/` - Contains Slurm job submission scripts at various HPC centers
3. `irfs/` - Contains IRF File lists for HPC centers that have a FEOTS database
3. `scripts/` - Contains post-processing scripts for analyzing simulation output

Underneath the `params/`,`slurm/`, and `irfs/` directories are subdirectories for namelist parameters, IRF file lists, and slurm job submission scripts specific to Google Cloud and LANL computing resources.

```
o zapiola
|
|
o --- o params/
|     |
|     |
o     o --- o google/
|     |
|     |
o     o --- o lanl/
|
|
o --- o slurm/
|     |
|     |
o     o --- o google/
|     |
|     |
o     o --- o lanl/
|
|
o --- o irfs/
|     |
|     |
o     o --- o google/
|     |
|     |
o     o --- o lanl/
|
|
o --- o scripts/
```

This example also contains
1. `FEOTSInitialize.f90` - A program for generating the initial conditions for the Zapiola simulation
2. `GenMask.f90` - A program for generating the regional mask for the Zapiola simulation
3. `makefile` - A makefile for building the `init` and `genmask` programs

## Running this demo

### On Turquoise
Set the following environment variables
```
export FEOTS_DBROOT=/lustre/scratch3/turquoise/jschoonover/feots/E3SMV0-HILAT-5DAVG
export OUTDIR=/lustre/scratch3/turquoise/jschoonover/zapiola
```
The `FEOTS_DBROOT` is the path to a FEOTS database. `OUTDIR` is the path to where FEOTS output is written.
1. Create the region mask - ` sbatch --get-user-env slurm/lanl/01-genmask.slurm`
2. Extract the regional operators from the global operators - ` sbatch --get-user-env slurm/lanl/02-regional-extraction.slurm`
3. Create the initial conditions - ` sbatch --get-user-env slurm/lanl/03-init.slurm`
4. Forward step the model - ` sbatch --get-user-env slurm/lanl/04-integrate.slurm`

To modify this experiment on Turquoise, you can edit `params/lanl/runtime.params`

### On a local workstation
Extraction of the regional operators requires about 30 GB RAM. If you have extracted the regional operators already, this regional simulation is lightweight enough that it can be run on a local workstation with 1 core/hyperthread and about 3 GB RAM.
Set the following environment variables
```
export REGIONAL_DB=/home/joe/feots_output/regional_constant_preservation
export OUTDIR=/home/joe/feots_output/regional_constant_preservation/no-mixing
```
If you would like to disable vertical mixing
```
export FEOTS_FLAGS='--no-vertical-mixing'
```
To unset the `FEOTS_FLAGS` variable, use `unset FEOTS_FLAGS`.

Make sure that the following files are present in the $OUTDIR
* `mappings.regional`
* `mask.nc`
* `mesh.nc`
* `transport.*.h5`
* `diffusion.*.h5`

If you plan on running the post-processor, make sure you install the python dependencies before running the pipeline

When you are ready, simply run `bash ./bash/job-pipeline.sh`. This will run through initial condition creation, forward integration, and post-processing stages. If you want to run the individual steps in the job pipeline, take the time to review the `bash/job-pipeline.sh` script to see how it works.


