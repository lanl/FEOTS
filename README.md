# FEOTs
Fast Equilibration of Ocean Tracers


## Getting Started

### Package Dependencies
* Fortran Compiler (2008 Compliant)
* MPI
* Parallel HDF5 ( C and Fortran builds )
* NetCDF-C
* NetCDF-Fortran

**Spack Dependency Install**
```
spack install netcdf-fortran ^netcdf-c+mpi ^hdf5+fortran
```

### Installation
```
autoreconf --install
FCFLAGS="-Ofast -ffree-line-length-none" ./configure --enable-mpi --prefix=/path/to/install
make
make install
```
