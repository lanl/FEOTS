lib_LIBRARIES = libfeots.a
bin_PROGRAMS = feots

AM_FCFLAGS = $(HDF5_FFLAGS) $(NETCDF4_CFLAGS) $(NETCDF4_LDFLAGS) -lnetcdff -lnetcdf
AM_LDFLAGS = $(HDF5_LDFLAGS) $(NETCDF4_LDFLAGS)
AM_CPPFLAGS = $(HDF5_CPPLAGS) $(NETCDF4_CPPFLAGS)

DISTCLEANFILES = modelprecision.mod constantsdictionary.mod commonroutines.mod binaryio.mod crsmatrix_class.mod tracerstorage_class.mod pop_adjacencygraph_class.mod pop_gridtypemappings.mod pop_params_class.mod pop_mesh_class.mod pop_native_class.mod pop_stencil_class.mod pop_regional_class.mod pop_feots_class.mod feots_driver_routines.mod feots.mod

include_HEADERS = modelprecision.mod constantsdictionary.mod commonroutines.mod binaryio.mod crsmatrix_class.mod tracerstorage_class.mod pop_adjacencygraph_class.mod pop_gridtypemappings.mod pop_params_class.mod pop_mesh_class.mod pop_native_class.mod pop_stencil_class.mod pop_regional_class.mod pop_feots_class.mod feots_driver_routines.mod

libfeots_a_SOURCES = ModelPrecision.F90 ConstantsDictionary.F90 CommonRoutines.F90 BinaryIO.F90 CRSMatrix_Class.F90 TracerStorage_Class.F90 POP_Stencil_Class.F90 POP_Mesh_Class.F90 POP_GridTypeMappings.F90 POP_AdjacencyGraph_Class.F90 POP_Params_Class.F90 POP_Native_Class.F90 POP_Regional_Class.F90 POP_FEOTS_Class.F90 FEOTS_Driver_Routines.F90

feots_SOURCES = ModelPrecision.F90 ConstantsDictionary.F90 CommonRoutines.F90 BinaryIO.F90 CRSMatrix_Class.F90 TracerStorage_Class.F90 POP_Stencil_Class.F90 POP_Mesh_Class.F90 POP_GridTypeMappings.F90 POP_AdjacencyGraph_Class.F90 POP_Params_Class.F90 POP_Native_Class.F90 POP_Regional_Class.F90 POP_FEOTS_Class.F90 FEOTS_Driver_Routines.F90 FEOTS.F90
