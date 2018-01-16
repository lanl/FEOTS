! GreedyGraphColoring.f90
!
! Author : Joseph Schoonover
! E-mail : jschoonover@lanl.gov, schoonover.numerics@gmail.com
!
! Copyright 2017 Joseph Schoonover, Los Alamos National Laboratory
! 
! Redistribution and use in source and binary forms, with or without
! modification,
! are permitted provided that the following conditions are met:
! 
! 1. Redistributions of source code must retain the above copyright notice, this
! list of conditions and the following disclaimer.
! 
! 2. Redistributions in binary form must reproduce the above copyright notice,
! this list of conditions and the following disclaimer in the 
! documentation and/or other materials provided with the distribution.
! 
! 3. Neither the name of the copyright holder nor the names of its contributors
! may be used to endorse or promote products derived from this 
!  software without specific prior written permission.
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED 
! TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
! PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR 
! CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
! EXEMPLARY,  OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
! BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
! LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE  OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! ////////////////////////////////////////////////////////////////////////////////////////////////

PROGRAM GreedyGraphColoring

! src/common/
!USE ModelPrecision
! src/POP/
USE POP_Params_Class
USE POP_Mesh_Class
USE POP_Stencil_Class
USE POP_AdjacencyGraph_Class
USE POP_Native_Class

IMPLICIT NONE


   TYPE( POP_Params )         :: params
   TYPE( POP_Mesh )           :: mesh
   TYPE( Stencil )            :: overlapStencil
   TYPE( POP_AdjacencyGraph ) :: graph
   TYPE( POP_Native )         :: impulseFields

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


      CALL impulseFields % Build( mesh, graph % nColors )
      CALL impulseFields % InitializeForNetCDFWrite( ImpulseField, &
                                                     mesh, &
                                                     'ImpulseFields.nc', &
                                                     .TRUE. )

      CALL GraphToImpulse( graph, impulseFields, mesh )

      CALL impulseFields % WriteTracerToNetCDF( mesh )
 
      CALL impulseFields % FinalizeNetCDF( )

      ! Write the graph to file for later use
      CALL graph % WriteGraphBinFile( TRIM(params % GraphFile) )


      ! Clear memory
      CALL overlapStencil % Trash( )
      CALL graph % Trash( )
      CALL mesh % Trash( )

CONTAINS

  SUBROUTINE GraphToImpulse( graph, impulse, mesh )
    IMPLICIT NONE
    TYPE( POP_AdjacencyGraph ), INTENT(in) :: graph
    TYPE( POP_Native ), INTENT(inout)      :: impulse
    TYPE( POP_Mesh ), INTENT(in)           :: mesh
    ! Local
    INTEGER :: irf_id, col, i, j, k


      DO irf_id = 1, graph % nColors
        DO col = 1, mesh % nDOF

          IF( graph % color(col) == irf_id )THEN

            i = mesh % dofToIJK(1,col)
            j = mesh % dofToIJK(2,col)
            k = mesh % dofToIJK(3,col)

            impulse % tracer(i,j,k,irf_id) = 1.0_prec

          ENDIF

        ENDDO
      ENDDO


  END SUBROUTINE GraphToImpulse

END PROGRAM GreedyGraphColoring

