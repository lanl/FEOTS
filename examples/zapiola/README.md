# FEOTS - Zapiola Demo

## Description


## Running this demo
Make sure you have the FEOTS toolkit installed. The install path is referred to as ${FEOTS_INSTALL}, and it is assumed that `$FEOTS_INSTALL/bin` is in your default search path.

Modify the `FEOTS_Settings` file to fit your system, then
```
$ module load <compiler> <mpi> <hdf5> <netcdf>
```
Upate the makefile in this directory to use the appropriate Fortran compiler (FC), and LIB and INCLUDE variables for NetCDF, HDF5, and FEOTS.

```
$ make all
```

To run this demo (assuming you have global IRFs),
```
$ ./GenMask
$ feots region-extraction
$ export OMP_NUM_THREADS=8
$ mpirun -np 7 -x OMP_NUM_THREADS ./init
$ mpirun -np 7 -x OMP_NUM_THREADS feots integrate
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



