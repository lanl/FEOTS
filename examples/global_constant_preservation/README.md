# FEOTS - Global Constant Preservation

This demo runs a simple global offline tracer simulation using the `E3SMV0-HILAT-5DAVG` FEOTS dataset. In this example, we aim to show that constant tracer fields are preserved when under the action of the transport (advection + lateral diffusion) and the vertical diffusion operators. The parent model (E3SMv0) uses finite volume discretizations of the advection and diffusion operators that are constant preserving. Thus, this example serves as a verification test that demonstrates FEOTS maintains certain properties of the parent model.


## Getting Started
In order to run this demonstration, you need the following packages installed on your system
* OpenMPI, MVAPICH, or MPICH
* HDF5 (Parallel build)
* NetCDF-C, NetCDF-Fortran
* FEOTS, configured with `--enable-mpi`

On LANL Turquoise network, we use `gcc/8.3.0`, `mvapich2/2.3`, `hdf5-parallel/1.8.16`, and `netcdf-h5parallel/4.4.0` to build FEOTS.

*Create the Initial Conditions*
To demonstrate constant preservation, we initialize a single tracer field with a constant value ( `tracer = 1` ) everywhere. To create these initial conditions, first compile the `FEOTSInitialize.f90` program provided.

```
$ module use /users/jschoonover/modulefiles
$ module load gcc/8.3.0 mvapich2/2.3 hdf5-parallel/1.8.16 netcdf-h5parallel/4.4.0 feots/dev
$ make init
```

Then, set the `FEOTS_DBROOT` and `OUTDIR` environment variables. 
`FEOTS_DBROOT` is the path to the FEOTS operator database. `OUTDIR` is the directory where the initial conditions will be written.

```
$ export FEOTS_DBROOT=/lustre/scratch3/turquoise/jschoonover/feots/E3SMV0-HILAT-5DAVG
$ export OUTDIR=/lustre/scratch3/turquoise/${USER}/global_constant_preservation
```

Use the `01-init.slurm` Slurm batch file to run the initial condition generator.
```
$ sbatch --get-user-env slurm/lanl/01-init.slurm
```
Notice that we pass the `--get-user-env` flag to `sbatch` so that your environment variables are passed to the environment on the compute nodes where the batch script is executed.

When this job is complete, your initial condition files can be found at `${OUTDIR}/Tracer.init.nc`.

*Forward Step the Model*
Now that you have initial conditions, you can run the model in a forward-integration mode. Forward-integration simply will forward-step the offline tracer model with the time stepping parameters indicated in `params/lanl/runtime.params`. 

Use the `02-integrate.slurm` Slurm batch file  to forward step the model.
```
$ sbatch --get-user-env slurm/lanl/02-integrate.slurm
```

When the job completes you will have the solution at one time-step later than the initial conditions. This output can be found at `${OUTDIR}/Tracer.0000000000.nc`

## Validate the results
You can confirm that constant preservation is achieved by comparing the initial condition file and the forward-stepped solution. We expect that the tracer, for all time, is constant; thus, the initial condition is the exact solution. Taking the difference of the forward-stepped solution with the initial condition then gives a measure of the error in the transport and diffusion in FEOTS.

Use the NCO routines to calculate the error,
```
$ ncdiff --variable Tracer_01 ${OUTDIR}/Tracer.0000000000.nc ${OUTDIR}/Tracer.init.nc ${OUTDIR}/Error.nc
```
Now, calculate the absolute maximum of the error using `ncwa`
```
$ ncwa --operation mabs --variable Tracer_01 ${OUTDIR}/Error.nc ${OUTDIR}/mabsError_T01.nc
```
Finally, you can use ncdump to print the maximum absolute value of the error.
```
$ ncdump --variable Tracer_01 ${OUTDIR}/mabsError_T01.nc
```
