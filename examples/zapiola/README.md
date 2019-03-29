# FEOTS - Zapiola Demo

## Description


## Running this demo

Modify the `FEOTS_Settings` file to fit your system, then
```
$ module load <compiler> <mpi> <hdf5> <netcdf>
$ source FEOTS_Settings
$ make all
```

To run this demo (assuming you have global IRFs),
```
$ ./GenMask
$ ./RegionalExtraction
$ export OMP_NUM_THREADS=8
$ mpirun -np 7 -x OMP_NUM_THREADS ./FEOTSInitialize
$ mpirun -np 7 -x OMP_NUM_THREADS ./FEOTSDriver
```

For your convenience, `feots.slurm` is a batch submission file for slurm for this demo.


## Modifying this demo
The `runtime.params` namelist file controls most of the settings for this demo, including
the lat/lon boundaries, time stepping settings, and location of your global and regional
IRF files.

### Mask Generation and Boundary Conditions
Boundary conditions for this regional simulation are controlled both by the GenMask.f90
program and the FEOTSInitialize.f90 program. The mask determines which lat/lon cells
are in the domain and are active, prescribed, or open. Active cells are forward stepped
in the driver. Prescribed cells are forced to retain the value initially prescribed to them.
Open cells are prescribed cells with a concentration value of 0.

The regional mask is generated in the GenMask.f90 program included in this sub-directory. The
prescribed cells are set to be the cells along the south, east, and north boundaries of this
domain. The width of the prescribed region is set in the `prescribed_width` parameter in the
GenMask.f90 program. Currently, it is set to 0.5 degrees. 

**If you change your mask at any time, you will have to re-run the RegionalExtraction program**



