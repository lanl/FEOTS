PROGRAM GreedyGraphColoring

! src/common/
!USE ModelPrecision
! src/POP/
USE POP_Params_Class
USE POP_Mesh_Class
USE POP_Stencil_Class
USE POP_AdjacencyGraph_Class

IMPLICIT NONE


   TYPE( POP_Params )         :: params
   TYPE( POP_Mesh )           :: mesh
   TYPE( Stencil )            :: overlapStencil
   TYPE( POP_AdjacencyGraph ) :: graph

      CALL params % Build( )

      ! Load in the mesh from the netcdf file specified above
      CALL mesh % Load( TRIM(params % meshfile) )
      mesh % meshType = params % meshType ! And set a flag for the type of mesh

      ! Here, we build a overlap stencil for the finite difference scheme
      ! specified above
      CALL overLapStencil % Build( stencilFlag = params % StencilType, &
                                   flavor      = Overlap )


      ! Build the adjacency graph using the wet points in the mesh and the
      ! overlap-stencil associated with the specified transport operator finite
      ! difference stencil
      CALL graph % ConstructFromStencil( mesh, overlapStencil ) 
      

      ! Now that we have an adjacency graph, we can now perform the greedy
      ! coloring. This coloring can be used to set initial tracer fields for
      ! diagnosing transport operators with the associated finite differencing
      ! scheme.
      CALL graph % GreedyColoring( )

      ! Write the graph to file for later use
      !CALL graph % WriteGraphFile( 'laxwendroff.pop.ptgrid.graph' )
      CALL graph % WriteGraphBinFile( TRIM(params % GraphFile) )


      ! Clear memory
      CALL overlapStencil % Trash( )
      CALL graph % Trash( )
      CALL mesh % Trash( )

END PROGRAM GreedyGraphColoring

