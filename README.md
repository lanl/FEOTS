# FEOTs
Fast Equilibration of Ocean Tracers

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5576912.svg)](https://doi.org/10.5281/zenodo.5576912)


## Getting Started
* [Run the Argentine Basin offline tracer simulation on Google Cloud](https://fluidnumerics.github.io/FEOTS/codelabs/feots-on-google-cloud/#0)

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
