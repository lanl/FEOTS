"""
HPC Base image

Contents:
  HDF5 version 1.10.5
  GNU compilers version 8.3
  OpenMPI version 4.0.1
  Python 2 and 3 (upstream)
"""
# pylint: disable=invalid-name, undefined-variable, used-before-assignment

devel_image = 'centos-7'
runtime_image = 'centos-7'

######
# Devel stage
######

Stage0 += comment(__doc__, reformat=False)

Stage0 += baseimage(image=devel_image, _as='devel')

# Python
Stage0 += python()

# GNU compilers
compiler = gnu(version='8')
Stage0 += compiler

# OpenMPI
Stage0 += openmpi(version='4.0.1', cuda=True, infiniband=False, toolchain=compiler.toolchain)

# HDF5
Stage0 += hdf5(version='1.10.5', mpi=True, toolchain=compiler.toolchain)

# NetCDF
Stage0 += netcdf(version='4.7.3', version_fortran='4.5.2', fortran=True, with_mpi=True, toolchain=compiler.toolchain)
