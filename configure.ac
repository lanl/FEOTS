AC_INIT([feots],[2.0.0],[support@fluidnumerics.com])
AM_INIT_AUTOMAKE(-Wall -Werror foreign)

# Set up build system to generate libSELF.a
AC_PROG_RANLIB
AM_PROG_AR

AC_CONFIG_MACRO_DIRS([m4])

# Find Fortran compilers
AC_PROG_FC([gfortran])

AC_LANG_PUSH([Fortran])

# Check for OpenMP
AC_ARG_ENABLE([openmp], AS_HELP_STRING([--enable-openmp], [default:no]))
AS_IF([test "x$enable_openmp" = "xyes"],[

  FCFLAGS+=" -DHAVE_OPENMP"
  AX_CHECK_COMPILE_FLAG([-mp], [FCFLAGS+=" -mp"])
  AX_CHECK_COMPILE_FLAG([-openmp], [FCFLAGS+=" -openmp"])
  AX_CHECK_COMPILE_FLAG([-fopenmp], [FCFLAGS+=" -fopenmp"])

])

# Check for HDF5
AX_LIB_HDF5()
AS_IF([test "x$with_hdf5" = "xyes"],[
  AC_MSG_NOTICE([found hdf5 $H5FC])
],[
  AC_MSG_ERROR([Unable to find HDF5.])
])

# Check for NetCDF4
AX_LIB_NETCDF4()
AS_IF([test "x$with_netcdf" = "xyes"],[
  FC="$NETCDF4_FC"
],[
  AC_MSG_ERROR([Unable to find NetCDF4.])
])

# Check for MPI
AC_CHECK_PROGS(MPIFC, mpif90 mpxlf95_r mpxlf90_r mpxlf95 mpxlf90 mpf90 cmpif90c, $FC)
AC_ARG_ENABLE([mpi], AS_HELP_STRING([--enable-mpi], [default:no]))
AS_IF([test "x$enable_mpi" = "xyes"],[
  FCFLAGS+=" -DHAVE_MPI"
  AS_IF([test "x$HDF5_TYPE" = "xserial",],[
    AC_MSG_ERROR([HDF5 found is serial. We need parallel HDF5.])
  ])

])

# Check for Double Precision floating point arithmetic request.
AC_ARG_ENABLE( double, AS_HELP_STRING([--enable-double, default:no]))
AS_IF([test "x$enable_double" = "xyes"],[
  FCFLAGS+=" -DDOUBLE_PRECISION"
])

AC_CONFIG_FILES([
Makefile
src/Makefile
])

AC_OUTPUT
