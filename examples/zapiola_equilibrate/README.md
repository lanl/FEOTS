# FEOTS - Zapiola Demo

## Description



## Directory Structure
This example contains directories for 
1. `slurm/` - Contains Slurm job submission scripts at various HPC centers
2. `irfs/` - Contains IRF File lists for HPC centers that have a FEOTS database

Underneath the `slurm/` and `irfs/` directories are subdirectories for namelist parameters, IRF file lists, and slurm job submission scripts specific to Google Cloud and LANL computing resources.

```
o zapiola
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
export REGIONAL_DB=/lustre/scratch3/turquoise/jschoonover/zapiola
export OUTDIR=/lustre/scratch3/turquoise/jschoonover/zapiola/simulation-1
```
The `FEOTS_DBROOT` is the path to a FEOTS database. `OUTDIR` is the path to where FEOTS output is written.
1. Create the region mask - ` sbatch --get-user-env slurm/lanl/01-genmask.slurm`
2. Extract the regional operators from the global operators - ` sbatch --get-user-env slurm/lanl/02-regional-extraction.slurm`
3. Create the mappings.regional file from the mask and the global operators - ` sbatch --get-user-env slurm/lanl/03-genmaps.slurm`
4. Create the initial conditions - ` sbatch --get-user-env slurm/lanl/04-init.slurm`
5. Forward step the model - ` sbatch --get-user-env slurm/lanl/05-integrate.slurm`

To modify this experiment on Turquoise, you can edit `runtime.params`
