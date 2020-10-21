! FEOTS_Driver_Routines.f90
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

MODULE FEOTS_Driver_Routines

USE ModelPrecision
USE CommonRoutines
USE CRSMatrix_Class
USE POP_Params_Class
USE POP_Mesh_Class
USE POP_Stencil_Class
USE POP_AdjacencyGraph_Class
USE POP_Native_Class
USE POP_FEOTS_Class
USE POP_GridTypeMappings
USE FEOTS_CLI_Class
#ifdef HAVE_OPENMP
USE OMP_LIB
#endif


IMPLICIT NONE

#include "FEOTS_Macros.h"


CONTAINS

  SUBROUTINE ExtractOceanState(cliParams)

    IMPLICIT NONE
    TYPE( FEOTS_CLI )    :: cliParams
    TYPE( POP_FEOTS )    :: feots
    TYPE( POP_Mesh )     :: globalMesh 
    TYPE( POP_Native )   :: globalState 
    
    INTEGER        :: fileID, fUnit
    INTEGER        :: m, i, j, k, i_local, j_local, k_local
    CHARACTER(5)   :: fileIDChar
    CHARACTER(200) :: thisIRFFile
    CHARACTER(200) :: oceanStateFile

     CALL feots % Build( cliParams, 0, 1 )

     CALL globalMesh % Load( TRIM(cliParams % dbRoot)//'/mesh/mesh.nc' )

     CALL globalState % Build( globalMesh, 1, 0, 1 ) 

     thisIRFFile = cliParams % irfFile
     PRINT*,' Loading '//TRIM(cliParams % irfFile)
     CALL globalState % LoadOceanState( globalMesh, cliParams % irfFile )
     ! The volume correction attribute is a unitless measure; it is
     ! the fractional change in the fluid volume at each degree of freedom
     ! The ssh divided by the volume's height is equivalent to the volume
     ! correction
     globalState % volume = globalState % volume/globalMesh % dzw(1)

     DO m = 1, feots % feotsMap % nCells
        i = feots % feotsMap % IJKinRegion(1,m)
        j = feots % feotsMap % IJKinRegion(2,m)
        k = feots % feotsMap % IJKinRegion(3,m)

        i_local = feots % feotsMap % dofToLocalIJK(1,m)
        j_local = feots % feotsMap % dofToLocalIJK(2,m)
        k_local = feots % feotsMap % dofToLocalIJK(3,m)

        feots % nativeSol % temperature(i_local,j_local,k_local) = globalState % temperature(i,j,k)
        feots % nativeSol % salinity(i_local,j_local,k_local)    = globalState % salinity(i,j,k)
        feots % nativeSol % density(i_local,j_local,k_local)     = globalState % density(i,j,k)
        !feots % nativeSol % volume(i_local,j_local,k_local)      = globalState % volume(i,j,k)
     ENDDO

     WRITE(fileIDChar, '(I5.5)' ) cliParams % oplevel
     oceanStateFile = TRIM( cliParams % outDir )//'Ocean.'//fileIDChar//'.nc'
     CALL feots % nativeSol % WriteOceanState( feots % mesh, TRIM(oceanStateFile) )

     CALL globalState % Trash( )
     CALL globalMesh % Trash( )
     CALL feots % Trash( )

  END SUBROUTINE ExtractOceanState

  SUBROUTINE GenerateMeshOnlyFile(cliParams)
   
   IMPLICIT NONE

   TYPE( FEOTS_CLI )  :: cliParams
   TYPE( POP_Params ) :: params
   TYPE( POP_Mesh )   :: mesh
   INTEGER(KIND=8), ALLOCATABLE :: dofToIJK_check(:,:), ijkToDOF_check(:,:,:)
   INTEGER :: i, j, k, fUnit, diffcount, thisdiff
   CHARACTER(400) :: ncfile


      ! This call loads the mesh from the netcdf file and builds the
      ! ijkToDOF and dofToIJK mappings by using the "kmt" field.
      CALL mesh % Load( TRIM(cliParams % irfFile)  )

      ! Write a netcdf file containing only the mesh
      CALL mesh % WriteNetCDF( TRIM(cliParams % dbRoot)//'/mesh/mesh.nc' )

      CALL mesh % Trash( )

  END SUBROUTINE GenerateMeshOnlyFile

  SUBROUTINE GreedyGraphColoring(cliParams)
    IMPLICIT NONE
    TYPE( FEOTS_CLI )          :: cliParams
    TYPE( POP_Params )         :: params
    TYPE( POP_Mesh )           :: mesh
    TYPE( Stencil )            :: overlapStencil
    TYPE( POP_AdjacencyGraph ) :: graph
    TYPE( POP_Native )         :: impulseFields

      CALL params % Build(cliParams % paramFile)

      ! Load in the mesh from the netcdf file specified above
      CALL mesh % Load( TRIM(cliParams % dbRoot)//'/mesh/mesh.nc')
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
      CALL graph % WriteGraph_HDF5( TRIM(cliParams %dbRoot)//'/irf/impulse/graph.h5' )

      CALL impulseFields % Build( mesh, graph % nColors, 0, 1 )
      CALL impulseFields % InitializeForNetCDFWrite( ImpulseField, &
                                                     mesh, &
                                                     TRIM(cliParams % dbRoot)//'/irf/impulse/ImpulseFields.nc', &
                                                     .TRUE. )

      CALL GraphToImpulse( graph, impulseFields, mesh )

      CALL impulseFields % WriteTracerToNetCDF( mesh )
 
      CALL impulseFields % FinalizeNetCDF( )

      ! Clear memory
      CALL impulseFields % Trash( )
      CALL overlapStencil % Trash( )
      CALL graph % Trash( )
      CALL mesh % Trash( )

  END SUBROUTINE GreedyGraphColoring

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

  SUBROUTINE GenMask(cliParams)

    IMPLICIT NONE
    TYPE( FEOTS_CLI )    :: cliParams
    TYPE( POP_Params )   :: params
    TYPE( POP_Mesh )     :: mesh
    INTEGER              :: i, j
    INTEGER, ALLOCATABLE :: maskField(:,:)
    CHARACTER(400)       :: ncfile
    REAL(prec)           :: x, y, r


      CALL params % Build(cliParams % paramFile)

      CALL mesh % Load( TRIM(cliParams % irfFile)  )

      ALLOCATE( maskField(1:mesh % nX,1:mesh % nY) )

      maskField = 0

      DO j = 1, mesh % nY
         DO i = 1, mesh % nX

            x = mesh % tLon(i,j)
            y = mesh % tLat(i,j)
          
            IF( x >= 180.0_prec )THEN
               x = x -360.0_prec
            ENDIF
            ! Build a circular region around the Agulhas            
            r = sqrt( (x-20.0_prec)**2 + (y+40.0_prec)**2 )
            IF( r <= 30.0_prec )THEN
               IF( r > 29.5_prec )THEN
                  maskfield(i,j) = -1 ! Prescribed Points
               ELSE
                  maskfield(i,j) = 1  ! Interior Points
               ENDIF
            ENDIF

         ENDDO
      ENDDO

      CALL WriteMaskField( mesh, maskField, TRIM(params % maskFile) )
      CALL mesh % Trash( )

 END SUBROUTINE GenMask

 SUBROUTINE WriteMaskField( mesh, maskfield, maskfile )

   IMPLICIT NONE
   TYPE( POP_Mesh ), INTENT(inout)      :: mesh
   INTEGER, INTENT(in)                  :: maskfield(1:mesh % nX, 1:mesh % nY)
   CHARACTER(*), INTENT(in)             :: maskfile
   ! Local
   INTEGER :: start(1:2), recCount(1:2)
   INTEGER :: ncid, varid, x_dimid, y_dimid

      start    = (/1, 1/)
      recCount = (/mesh % nX, mesh % nY/)

      CALL Check( nf90_create( PATH=TRIM(maskfile),&
                               CMODE=OR(nf90_clobber,nf90_64bit_offset),&
                               NCID=ncid ) )
      CALL Check( nf90_def_dim( ncid, "nlon", mesh % nX, x_dimid ) )
      CALL Check( nf90_def_dim( ncid, "nlat", mesh % nY, y_dimid ) )
      CALL Check( nf90_def_var( ncid, "mask",NF90_INT,&
                               (/ x_dimid, y_dimid /),&
                                varid ) )

      CALL Check( nf90_put_att( ncid, varid, "long_name", "Domain Mask" ) )
      CALL Check( nf90_put_att( ncid, varid, "units", "" ) )
!      CALL Check( nf90_put_att( ncid, varid, "_FillValue", fillValue) )
!      CALL Check( nf90_put_att( ncid, varid, "missing_value", fillValue) )

      CALL Check( nf90_enddef(ncid) )

      CALL Check( nf90_put_var( ncid, &
                                varid, &
                                maskfield, &
                                start, recCount ) )


      CALL Check( nf90_close( ncid ) )


 END SUBROUTINE WriteMaskField

 SUBROUTINE OperatorDiagnosis(cliParams)

   
   IMPLICIT NONE
   TYPE( FEOTS_CLI )          :: cliParams

   TYPE( POP_Params )         :: params
   TYPE( POP_Mesh )           :: mesh
   TYPE( Stencil )            :: advstencil
   TYPE( POP_AdjacencyGraph ) :: graph
   TYPE( CRSMatrix )          :: transportOp, diffusionOp
   TYPE( POP_IRF )            :: irfFields
   CHARACTER(200)             :: thisIRFFile, meshFile
   CHARACTER(200)             :: crsFile
   CHARACTER(17)              :: walltime
   CHARACTER(5)               :: fileIDChar
   INTEGER                    :: nEntries, irf_id, row, col, m, diff_id
   INTEGER                    :: i, j, k, this_i, this_j, this_k, iel
   INTEGER                    :: true_i, true_j, kdiff
   INTEGER                    :: fUnit, fileID
   INTEGER, ALLOCATABLE       :: nval(:)
   REAL(prec), ALLOCATABLE    :: x(:), Dx(:)
   REAL(prec)                 :: t0, t1    


      CALL params % Build(cliParams % paramFile)

      ! Load in the mesh from the netcdf file specified above
      CALL mesh % Load( TRIM(cliParams % dbRoot)//'/mesh/mesh.nc')
      mesh % meshType = params % meshType ! And set a flag for the type of mesh

      ! Here, we build a stencil for the finite difference scheme
      ! specified above, without overlaps
      CALL advStencil % Build( stencilFlag = params % StencilType, &
                               flavor      = Normal )

      ! Allocate space for the transport operator
      ! The number of possible non-zero entries is the number of degrees of
      ! freedom multiplied by the stencil-width
      nEntries = ( mesh % nDOF )*( advStencil % nPoints )
      CALL transportOp % Build( INT(mesh % nDOF,8), INT(mesh % nDOF,8), INT(nEntries,8) )
      CALL diffusionOp % Build( INT(mesh % nDOF,8), INT(mesh % nDOF,8), INT(mesh % nDOF*3,8) ) 

      ! Set the row boundary indices assuming that the number of columns per row
      ! is fixed
      CALL transportOp % SetRowBounds( maxColPerRow = advStencil % nPoints )
      CALL diffusionOp % SetRowBounds( maxColPerRow = 3 )

      ! Allocate space for a hash-table storage
      ALLOCATE( nval(1:mesh % nDOF) )
      ALLOCATE( x(1:mesh % nDOF), Dx(1:mesh % nDOF) )


      ! Read the adjacency graph from file
      CALL graph % ReadGraph_HDF5( TRIM(cliParams % dbRoot)//'/irf/impulse/graph.h5' )
      
      ! Allocate space for the IRF fields. The number of IRF fields is
      ! equivalent to the number of colors in the graph
      !
      ! When building this data structure for the irf-fields, nColors
      ! corresponds to the number of irf fields-we add an additionation
      ! "POP_Native" attribute for the vertical diffusivity
      CALL irfFields % Build( mesh % nX, &
                              mesh % nY, &
                              mesh % nZ, &
                              graph % nColors )
      
      ! /////////////////////// Op Diagnosis //////////////////////////// !

      fileID = cliParams % oplevel
      thisIRFFile = cliParams % irfFile
 
      CALL irfFields % LoadNetCDF( filename=TRIM(thisIRFFile) )
  
      PRINT*, 'Diagnosing Transport Ops.'
      CALL CPU_TIME( t0 )
      nval    = 0
      ! Each degree of freedom in the IRF gives a column or the transport
      ! matrix.
      DO col = 1, mesh % nDOF 

         ! Obtain the (i,j,k)-tuple for this IRF degree of freedom
         i = mesh % DOFToIJK(1,col)
         j = mesh % DOFToIJK(2,col)
         k = mesh % DOFToIJK(3,col)

         DO m = 1, advstencil % nPoints

            this_i = i + advstencil % relativeNeighbors(1,m)
            this_j = j + advstencil % relativeNeighbors(2,m)
            this_k = k + advstencil % relativeNeighbors(3,m)
         
            ! Get neighboring i,j indices, taking tripole grid connectivity into account
            CALL GetTrueIJ( params % MeshType, &
                            this_i, this_j, &
                            mesh % nX, mesh % nY, &
                            true_i, true_j )

            IF( this_k <= mesh % KMT(true_i, true_j) .AND. this_k > 0 )THEN

               row = mesh % IJKToDOF(true_i,true_j,this_k)

               IF( row /= 0 )THEN

                  irf_id = graph % color(col)
                  nval(row) = nval(row) + 1
                  iel = transportOp % rowBounds(1,row) + nval(row)-1
                  transportOp % A(iel) = irfFields % irf(true_i,true_j,this_k,irf_id) 
                  transportOp % col(iel) = col
                  
               ENDIF

            ENDIF

         ENDDO

      ENDDO
      CALL transportOp % CountZeroValues( )

      PRINT*, 'Diagnosing Vertical Diffusion Op.'

      diff_id = graph % nColors + 1
      ! Diffusion operator diagnosis
      DO row = 1, mesh % nDOF

         i = mesh % DOFToIJK(1,row)
         j = mesh % DOFToIJK(2,row)
         k = mesh % DOFToIJK(3,row)

         IF( mesh % KMT(i,j) > 1 .AND. k <= mesh % KMT(i,j) )THEN

           IF( k == 1 )THEN ! At the surface, the "no-diffusive-flux" condition is used
               ! cell k coefficient/diagonal
               iel = diffusionOp % rowBounds(1,row)
               diffusionOp % A(iel) = -irfFields % irf(i,j,k,diff_id)/mesh % dz(k)/mesh % dzw(k)
               diffusionOp % col(iel) = row

               ! cell k+1
               iel = diffusionOp % rowBounds(1,row) + 1
               diffusionOp % A(iel) = irfFields % irf(i,j,k,diff_id)/mesh % dz(k)/mesh % dzw(k)
               ! our DOF varies with k varying fastest
               diffusionOp % col(iel) = row+1 !mesh % IJKToDOF(i,j,k+1)
               IF( row == 1 )THEN
                 PRINT*, 'diffOp Diag : ', diffusionOp % rowBounds(1,row), &
                                           i,j,k, &
                                           irfFields % irf(i,j,k,diff_id), &
                                           diffusionOp % A(iel)
               ENDIF

           ELSEIF( k == mesh % KMT(i,j) )THEN !At the bottom-- no diffusive flux

               ! cell k-1
               iel = diffusionOp % rowBounds(1,row)
               diffusionOp % A(iel) = irfFields % irf(i,j,k-1,diff_id)/mesh % dz(k)/mesh % dzw(k-1)
               diffusionOp % col(iel) = row-1!mesh % IJKToDOF(i,j,k-1)

               ! cell k coefficient/diagonal
               iel = diffusionOp % rowBounds(1,row) + 1
               diffusionOp % A(iel) = -irfFields % irf(i,j,k-1,diff_id)/mesh % dz(k)/mesh % dzw(k-1)
               diffusionOp % col(iel) = row

           ELSE

               ! cell k-1
               iel = diffusionOp % rowBounds(1,row)
               diffusionOp % A(iel) = irfFields % irf(i,j,k-1,diff_id)/mesh % dz(k)/mesh % dzw(k-1)
               diffusionOp % col(iel) = row-1!mesh % IJKToDOF(i,j,k-1)

               ! cell k coefficient/diagonal
               iel = diffusionOp % rowBounds(1,row) + 1
               diffusionOp % A(iel) = -irfFields % irf(i,j,k-1,diff_id)/mesh % dz(k)/mesh % dzw(k-1)-&
                                                   irfFields % irf(i,j,k,diff_id)/mesh % dz(k)/mesh % dzw(k)
               diffusionOp % col(iel) = row

               ! cell k+1
               iel = diffusionOp % rowBounds(1,row) + 2
               diffusionOp % A(iel) = irfFields % irf(i,j,k,diff_id)/mesh % dz(k)/mesh % dzw(k)
               diffusionOp % col(iel) = row+1!mesh % IJKToDOF(i,j,k+1)

           ENDIF
         ENDIF                 
      ENDDO
      CALL CPU_TIME( t1 )
      PRINT*, 'DONE! '
      WRITE( walltime, '(F17.4)' ) t1-t0
      PRINT*, 'Time to diagnose operator : '//TRIM(walltime)//' sec'

      WRITE(fileIDChar, '(I5.5)' ) fileID
      crsFile=TRIM(cliParams % dbRoot)//'/ops/transport.'//fileIDChar//'.h5'
      PRINT*,'Writing CRS Matrix files : '//TRIM(crsFile)
      CALL transportOp % WriteCRSMatrix_HDF5( TRIM(crsFile) )

      crsFile=TRIM(cliParams % dbRoot)//'/ops/diffusion.'//fileIDChar//'.h5'
      PRINT*,'Writing CRS Matrix files : '//TRIM(crsFile)
      CALL diffusionOp % WriteCRSMatrix_HDF5( TRIM(crsFile) )

      ! Testing the diffusion operator
      PRINT*, 'Testing diffusion operator'
      x = 1.0_prec
      Dx = diffusionOp % MatVecMul( x )
      PRINT*, 'MAX RowSum(Dx)', MAXVAL(Dx)
      PRINT*, 'MIN RowSum(Dx)', MINVAL(Dx)

      !CALL transportOp % Reset( )
      !CALL diffusionOp % Reset( )

      !ENDIF

      !ENDDO ! Loop over the IRF Files
      !CLOSE( fUnit )

      ! Clear memory
      CALL advStencil % Trash( )
      CALL graph % Trash( )
      CALL mesh % Trash( )
      CALL irfFields % Trash( )
      CALL transportOp % Trash( )
      CALL diffusionOp % Trash( )
      DEALLOCATE( x, Dx )
      DEALLOCATE( nval )!, columns )
      !DEALLOCATE( opdata )

  END SUBROUTINE OperatorDiagnosis

  SUBROUTINE RegionalMaps(cliParams)

    IMPLICIT NONE
    TYPE( FEOTS_CLI )    :: cliParams
    TYPE( CRSMatrix )    :: transportOp, diffusionOP
    TYPE( CRSMatrix )    :: regionalTransportOp, regionalDiffusionOP
    TYPE( POP_Params )   :: params
    TYPE( POP_Mesh )     :: globalMesh 
    TYPE( POP_Mesh )     :: regionalMesh 
    TYPE( Stencil )      :: modelstencil, advStencil
    TYPE( POP_Regional ) :: region
    
    INTEGER        :: fileID
    INTEGER     ::  nGentries, nRentries
    !INTEGER, ALLOCATABLE :: maskfield(:,:)
    CHARACTER(5)   :: fileIDChar
    CHARACTER(400) :: crsFile
    INTEGER :: irfStart, irfEnd

      CALL params % Build(cliParams % paramFile)

      CALL globalMesh % Load( TRIM(cliParams % dbRoot)//'/mesh/mesh.nc'  )

      CALL modelstencil % Build( stencilFlag = params % stencilType, &
                                 flavor      = LateralPlusCorners )

      CALL region % Build( globalMesh, regionalMesh, modelstencil, params % meshType, &
                           params % south, params % north, &
                           params % east, params % west, &
                           TRIM(cliParams % outdir)//'/mask.nc' ) 

      CALL regionalMesh % WriteNetCDF( TRIM(cliParams % outdir)//'/mesh.nc' )
      ! Write the regional data structure to a pickup file for later use

      CALL region % WritePickup( TRIM(cliParams % outdir)//'/mappings', maskProvided=.TRUE. )

      ! Clean up memory !
      CALL globalMesh % Trash( )
      CALL regionalMesh % Trash( )
      CALL region % Trash( )    

  END SUBROUTINE RegionalMaps

  SUBROUTINE RegionalExtraction(cliParams)

    IMPLICIT NONE
    TYPE( FEOTS_CLI )    :: cliParams
    TYPE( CRSMatrix )    :: transportOp, diffusionOP
    TYPE( CRSMatrix )    :: regionalTransportOp, regionalDiffusionOP
    TYPE( POP_Params )   :: params
    TYPE( POP_Mesh )     :: globalMesh 
    TYPE( POP_Mesh )     :: regionalMesh 
    TYPE( Stencil )      :: modelstencil, advStencil
    TYPE( POP_Regional ) :: region
    
    INTEGER        :: fileID
    INTEGER     ::  nGentries, nRentries
    !INTEGER, ALLOCATABLE :: maskfield(:,:)
    CHARACTER(5)   :: fileIDChar
    CHARACTER(400) :: crsFile
    INTEGER :: irfStart, irfEnd

      CALL params % Build(cliParams % paramFile)

      CALL globalMesh % Load( TRIM(cliParams % dbRoot)//'/mesh/mesh.nc'  )

      CALL modelstencil % Build( stencilFlag = params % stencilType, &
                                 flavor      = LateralPlusCorners )

      CALL region % Build( globalMesh, regionalMesh, modelstencil, params % meshType, &
                           params % south, params % north, &
                           params % east, params % west, &
                           TRIM(cliParams % regionalDb)//'/mask.nc' ) 

      CALL regionalMesh % WriteNetCDF( TRIM(cliParams % regionalDb)//'/mesh.nc' )

      CALL advstencil % Build( stencilFlag = params % stencilType, &
                               flavor      = Normal )

      nGentries = ( globalMesh % nDOF )*( advStencil % nPoints )
      CALL transportOp % Build( INT(globalmesh % nDOF,8), INT(globalmesh % nDOF,8), INT(nGEntries,8) )
      CALL diffusionOp % Build( INT(globalmesh % nDOF,8), INT(globalmesh % nDOF,8), INT(globalmesh % nDOF*3,8) ) 
      
      nRentries = ( region % nCells )*( advStencil % nPoints )
      CALL regionalTransportOp % Build( INT(region % nCells,8), INT(region % nCells,8), INT(nREntries,8) )
      CALL regionalDiffusionOp % Build( INT(region % nCells,8), INT(region % nCells,8), INT(region % nCells*3,8) ) 

      PRINT*, ' Extracting regional operators.'

      ! offset the file-id by the oplevel
      WRITE(fileIDChar, '(I5.5)' ) cliParams % oplevel
      crsFile=TRIM(cliParams % dbRoot)//'/ops/transport.'//fileIDChar//'.h5'
      PRINT*,'Reading CRS Matrix files : '//TRIM(crsFile)
      CALL transportOp % ReadCRSMatrix_HDF5(TRIM(crsFile), 0, 1 )

      crsFile=TRIM(cliParams % dbRoot)//'/ops/diffusion.'//fileIDChar//'.h5'
      PRINT*,'Reading CRS Matrix files : '//TRIM(crsFile)
      CALL diffusionOp % ReadCRSMatrix_HDF5(TRIM(crsFile), 0, 1 )

      ! Extract advection operator in region 
      CALL transportOp % SubSample( regionalTransportOp, &
                                    region % dofInRegion, &
                                    region % inverseDOFMap, &
                                    region % nCells, &
                                    nRentries )
      ! Extract diffusion operator in region 
      CALL diffusionOp % SubSample( regionalDiffusionOp, &
                                    region % dofInRegion, &
                                    region % inverseDOFMap, &
                                    region % nCells, &
                                    region % nCells*3 )

      crsFile=TRIM(cliParams % regionalDb)//'/transport.'//fileIDChar//'.h5'
      PRINT*,'Writing CRS Matrix files : '//TRIM(crsFile)
      CALL regionalTransportOp % WriteCRSMatrix_HDF5( TRIM(crsFile) )

      crsFile=TRIM(cliParams % regionalDb)//'/diffusion.'//fileIDChar//'.h5'
      PRINT*,'Writing CRS Matrix files : '//TRIM(crsFile)
      CALL regionalDiffusionOp % WriteCRSMatrix_HDF5( TRIM(crsFile) )
         
      CALL transportOp % Trash( )
      CALL diffusionOp % Trash( )
      CALL regionaltransportOp % Trash( )
      CALL regionalDiffusionOp % Trash( )
        
      ! Clean up memory !
      CALL globalMesh % Trash( )
      CALL regionalMesh % Trash( )
      CALL modelstencil % Trash( )
      CALL region % Trash( )    

  END SUBROUTINE RegionalExtraction

  SUBROUTINE FEOTSInitialize(cliParams)
    
    IMPLICIT NONE
    TYPE( FEOTS_CLI ) :: cliParams
    TYPE( POP_FEOTS ) :: feots
    CHARACTER(200)    :: thisIRFFile
    INTEGER           :: fUnit
    INTEGER :: myRank, nProcs, mpiErr

#ifdef HAVE_MPI
      CALL MPI_INIT( mpiErr )
      CALL MPI_COMM_RANK( MPI_COMM_WORLD, myRank, mpiErr )
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD, nProcs, mpiErr )
#else
      myRank = 0
      nProcs = 1
#endif
      CALL feots % Build( cliParams, myRank, nProcs )

      CALL InitialConditions( feots )

      !  //////////////////////////////////////////// File I/O  //////////////////////////////////////////////////////// !
      CALL feots % nativeSol % InitializeForNetCDFWrite( feots % params % TracerModel, &
                                                         feots % mesh, &
                                                         TRIM(cliParams % outDir)//'/Tracer.init.nc', &
                                                         .TRUE. )
      CALL feots % nativeSol % WriteNetCDFRecord( feots % mesh, 1 )
      CALL feots % nativeSol % WriteSourceEtcNetCDF( feots % mesh )

      CALL feots % nativeSol % FinalizeNetCDF( )
      ! //////////////////////////////////////////////////////////////////////////////////////////////////////////////// !

      CALL feots % Trash( )
#ifdef HAVE_MPI
      CALL MPI_FINALIZE( mpiErr )
#endif

 END SUBROUTINE FEOTSInitialize

 SUBROUTINE InitialConditions( myfeots )
 ! Sets the initial tracer distributions
   IMPLICIT NONE
   TYPE( POP_FEOTS ), INTENT(inout) :: myFeots
   ! Local
   INTEGER  :: i, j, k
   REAL(prec) :: x, y, z


      DO k = 1, myFeots % mesh % nZ  

         z = myFeots % mesh % z(k)

         DO j = 1, myFeots % mesh % nY
            DO i = 1, myFeots % mesh % nX 
 
               x = myFeots % mesh % tLon(i,j)
               y = myFeots % mesh % tLat(i,j)

               myFeots % nativeSol % tracer(i,j,k,:)  = 1.0_prec
               myFeots % nativeSol % source(i,j,k,:)  = 0.0_prec
               myFeots % nativeSol % rFac(i,j,k,:)    = 0.0_prec
               myFeots % nativeSol % mask(i,j,k,:)    = 1.0_prec


            ENDDO
         ENDDO
      ENDDO

 END SUBROUTINE InitialConditions

 SUBROUTINE FEOTSIntegrate(cliParams)
#undef __FUNC__
#define __FUNC__ "FeotsIntegrate"
   IMPLICIT NONE
   TYPE( FEOTS_CLI ), INTENT(in) :: cliParams
   ! Local
   TYPE( POP_FEOTS ) :: feots
   CHARACTER(10) :: ncFileTag
   CHARACTER(5)  :: fileIDChar
   CHARACTER(5)  :: rankChar
   CHARACTER(200):: thisIRFFile
   INTEGER       :: funit, recordID, fileID, i, nIODumps
   INTEGER       :: mpiErr, myRank, nProcs, iter
   REAL(prec)    :: tn
   REAL(prec)    :: t1, t2
   
      CALL MPI_INIT( mpiErr )
      CALL MPI_COMM_RANK( MPI_COMM_WORLD, myRank, mpiErr )
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD, nProcs, mpiErr )

      CALL feots % Build( cliParams, myRank, nProcs )
      WRITE( rankChar, '(I5.5)' ) myRank

      recordID = 1

      fileID   = feots % params % iterInit
      ! /////////////////////////// Load in the initial conditions //////////////////// !
      WRITE( ncfileTag, '(I10.10)' ) fileID
      ! For now, the record ID is 1. In the future, this will need to be
      ! calculated as a function of the initial iterate, the dump frequency,
      ! and the number of records per netcdf file
      
      !! Tracer.init.nc is read for the mask and source terms
      CALL feots % nativeSol % InitializeForNetCDFRead( feots % params % TracerModel,TRIM(cliParams % outDir)//'/Tracer.'//rankChar//'.init.nc', .TRUE. )
      CALL feots % nativeSol % ReadSourceEtcNetCDF( feots % mesh )
      CALL feots % nativeSol % ReadNetCDFRecord( feots % mesh, recordID )
      CALL feots % nativeSol % FinalizeNetCDF( )


         !***********

         IF( feots % params % iterInit == 0 )THEN
            tn = 0.0_prec
         ELSE
            ! This pickup file is read for the correct "initial condition"
            CALL feots % nativeSol % InitializeForNetCDFRead( feots % params % TracerModel,TRIM(cliParams % outDir)//'/Tracer.'//rankChar//'.'//ncFileTag//'.nc', .FALSE. )
            CALL feots % nativeSol % ReadNetCDFRecord( feots % mesh, recordID )
            CALL feots % nativeSol % FinalizeNetCDF( )
            tn = REAL(feots % params % iterInit,prec)*feots % params % dt
         ENDIF
   
         ! /////////////////////////////////////////////////////////////////////////////// !
         ! Transfer the data from the native storage to the FEOTS storage
         CALL feots % MapAllToDOF( )
         ! //// Forward Mode //// !
         PRINT*, '  Starting ForwardStep'
   
         DO iter = feots % params % iterInit, feots % params % iterInit + feots % params % nTimeSteps -1, feots % params % nStepsPerDump

#ifdef _OPENMP
            t1 = omp_get_wtime()
#else
            CALL CPU_TIME(t1)
#endif
            CALL feots % ForwardStep( tn, feots % params % nStepsPerDump, myRank, nProcs )
#ifdef _OPENMP
            t2 = omp_get_wtime()
#else
            CALL CPU_TIME(t2)
#endif
            INFO('ForwardStep wall-time = '//Float2Str(t2-t1)//' s')
 
            CALL feots % MapTracerFromDOF( )

            WRITE( ncfileTag, '(I10.10)' ) iter + feots % params % nStepsPerDump
            CALL feots % nativeSol % InitializeForNetCDFWrite( feots % params % TracerModel, &
                                                               feots % mesh, &
                                                               TRIM(cliParams % outdir)//'/Tracer.'//rankChar//'.'//ncFileTag//'.nc', &
                                                               .FALSE. )
            CALL feots % nativeSol % WriteNetCDFRecord( feots % mesh, recordID )
            CALL feots % nativeSol % FinalizeNetCDF( )

         ENDDO

      !CALL MPI_BARRIER( MPI_COMM_WORLD )
      CALL feots % Trash( )

      CALL MPI_FINALIZE( mpiErr )

 END SUBROUTINE FEOTSIntegrate

 SUBROUTINE FEOTSEquilibrate(cliParams)

   IMPLICIT NONE
   TYPE( FEOTS_CLI ) :: cliParams
   TYPE( POP_FEOTS ) :: feots
   CHARACTER(10) :: ncFileTag
   CHARACTER(5)  :: fileIDChar
   CHARACTER(200):: thisIRFFile
   INTEGER       :: funit, recordID, fileID, i, nIODumps
   INTEGER       :: mpiErr, myRank, nProcs, iter
   REAL(prec)    :: tn
   REAL(prec)    :: t1, t2
   
#ifdef HAVE_MPI
      CALL MPI_INIT( mpiErr )
      CALL MPI_COMM_RANK( MPI_COMM_WORLD, myRank, mpiErr )
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD, nProcs, mpiErr )
#else
      myRank = 0
      nProcs = 1
#endif

      CALL feots % Build( cliParams, myRank, nProcs )

      recordID = 1
      IF( myRank == 0 )THEN
        ! Tracer.init.nc is read for the mask and source terms
         CALL feots % nativeSol % InitializeForNetCDFRead( feots % params % TracerModel,'Tracer.init.nc', .TRUE. )
         CALL feots % nativeSol % ReadSourceEtcNetCDF( feots % mesh )
         CALL feots % nativeSol % ReadNetCDFRecord( feots % mesh, recordID )
         CALL feots % nativeSol % FinalizeNetCDF( )
   
  
         IF( feots % params % isPickupRun )THEN  
            CALL feots % nativeSol % InitializeForNetCDFRead( feots % params % TracerModel,'Tracer.pickup.nc', .FALSE. )
            CALL feots % nativeSol % ReadNetCDFRecord( feots % mesh, recordID )
            CALL feots % nativeSol % FinalizeNetCDF( )
            tn = REAL(feots % params % iterInit,prec)*feots % params % dt
   
         ENDIF
   
         ! /////////////////////////////////////////////////////////////////////////////// !
         
         ! Transfer the data from the native storage to the FEOTS storage
         CALL feots % MapAllToDOF( )
      ENDIF

#ifdef HAVE_OPENMP   
      t1 = omp_get_wtime( )
#else
      CALL CPU_TIME( t1 )
#endif
      CALL feots % JFNK( myRank )

          
#ifdef HAVE_OPENMP
      t2 = omp_get_wtime( )
#else
      CALL CPU_TIME( t2 )
#endif
      PRINT*, 'JFNK wall time :', t2-t1

      CALL feots % Trash( )
#ifdef HAVE_MPI
      CALL MPI_FINALIZE( mpiErr )
#endif

  END SUBROUTINE FEOTSEquilibrate

END MODULE FEOTS_Driver_Routines
