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
#ifdef HAVE_OPENMP
USE OMP_LIB
#endif


IMPLICIT NONE
CONTAINS
        
  SUBROUTINE ExtractOceanState()


    IMPLICIT NONE

    TYPE( POP_FEOTS )    :: feots
    TYPE( POP_Mesh )     :: globalMesh 
    TYPE( POP_Native )   :: globalState 
    
    INTEGER        :: fileID, fUnit
    INTEGER        :: m, i, j, k, i_local, j_local, k_local
    CHARACTER(5)   :: fileIDChar
    CHARACTER(200) :: thisIRFFile
    CHARACTER(200) :: oceanStateFile

      CALL feots % Build( 0, 1 )

      CALL globalMesh % Load( TRIM( feots % params % meshFile ) )

      CALL globalState % Build( globalMesh, 1 ) 

      OPEN( UNIT=NewUnit(fUnit),&
            FILE=TRIM(feots % params % IRFListFile), &
            FORM='FORMATTED',&
            ACCESS='SEQUENTIAL',&
            ACTION='READ',&
            STATUS='OLD' )

      DO fileID = 1, feots % params % nIRFFiles

         READ( fUnit, '(A200)' ) thisIRFFile

         IF( fileID >= feots % params % IRFStart )THEN
            PRINT*,' Loading '//TRIM(thisIRFFile)
            CALL globalState % LoadOceanState( globalMesh, thisIRFFile )
            ! The volume correction attribute is a unitless measure; it is
            ! the fractional change in the fluid volume at each degree of freedom
            ! The ssh divided by the volume's height is equivalent to the volume
            ! correction
            globalState % volume = globalState % volume/globalMesh % dzw(1)

            DO m = 1, feots % regionalMaps % nCells
               i = feots % regionalMaps % IJKinRegion(1,m)
               j = feots % regionalMaps % IJKinRegion(2,m)
               k = feots % regionalMaps % IJKinRegion(3,m)

               i_local = feots % regionalMaps % dofToLocalIJK(1,m)
               j_local = feots % regionalMaps % dofToLocalIJK(2,m)
               k_local = feots % regionalMaps % dofToLocalIJK(3,m)

               feots % nativeSol % temperature(i_local,j_local,k_local) = globalState % temperature(i,j,k)
               feots % nativeSol % salinity(i_local,j_local,k_local)    = globalState % salinity(i,j,k)
               feots % nativeSol % density(i_local,j_local,k_local)     = globalState % density(i,j,k)
               feots % nativeSol % volume(i_local,j_local,k_local)      = globalState % volume(i,j,k)
            ENDDO

            WRITE(fileIDChar, '(I5.5)' ) fileID
            oceanStateFile = TRIM( feots % params % regionalOperatorDirectory )//'Ocean.'//fileIDChar//'.nc'
            CALL feots % nativeSol % WriteOceanState( feots % mesh, TRIM(oceanStateFile) )

         ENDIF

     ENDDO

     CALL globalState % Trash( )
     CALL globalMesh % Trash( )
     CALL feots % Trash( )

  END SUBROUTINE ExtractOceanState

  SUBROUTINE GenerateMeshOnlyFile()
   
   IMPLICIT NONE

   TYPE( POP_Params ) :: params
   TYPE( POP_Mesh )   :: mesh
   INTEGER(KIND=8), ALLOCATABLE :: dofToIJK_check(:,:), ijkToDOF_check(:,:,:)
   INTEGER :: i, j, k, fUnit, diffcount, thisdiff
   CHARACTER(400) :: ncfile


      CALL params % Build( )

      ! Reads in the first file from the IRF File list
      OPEN( UNIT=NewUnit(fUnit),&
            FILE=TRIM(params % IRFListFile), &
            FORM='FORMATTED',&
            ACCESS='SEQUENTIAL',&
            ACTION='READ',&
            STATUS='OLD' )
      READ( fUnit, '(A400)' ) ncFile
      CLOSE( fUnit )

      ! This call loads the mesh from the netcdf file and builds the
      ! ijkToDOF and dofToIJK mappings by using the "kmt" field.
      CALL mesh % Load( TRIM(ncFile)  )

      ! Write a netcdf file containing only the mesh
      CALL mesh % WriteNetCDF( TRIM(params % meshFile) )


      CALL mesh % Trash( )

  END SUBROUTINE GenerateMeshOnlyFile

  SUBROUTINE GreedyGraphColoring()


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


    !  CALL impulseFields % Build( mesh, graph % nColors )
    !  CALL impulseFields % InitializeForNetCDFWrite( ImpulseField, &
    !                                                 mesh, &
    !                                                 'ImpulseFields.nc', &
    !                                                 .TRUE. )

    !  CALL GraphToImpulse( graph, impulseFields, mesh )

    !  CALL impulseFields % WriteTracerToNetCDF( mesh )
 
    !  CALL impulseFields % FinalizeNetCDF( )

      ! Write the graph to file for later use
      CALL graph % WriteGraphBinFile( TRIM(params % GraphFile) )


      ! Clear memory
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

  SUBROUTINE GenMask()

    IMPLICIT NONE

    TYPE( POP_Params )   :: params
    TYPE( POP_Mesh )     :: mesh
    INTEGER              :: i, j
    INTEGER, ALLOCATABLE :: maskField(:,:)
    CHARACTER(400)       :: ncfile
    REAL(prec)           :: x, y, r


      CALL params % Build( )

      CALL mesh % Load( TRIM(params % meshFile)  )

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

 SUBROUTINE OperatorDiagnosis()

   
   IMPLICIT NONE

   TYPE( POP_Params )         :: params
   TYPE( POP_Mesh )           :: mesh
   TYPE( Stencil )            :: advstencil
   TYPE( POP_AdjacencyGraph ) :: graph
   TYPE( CRSMatrix )          :: transportOperator, diffusionOperator
   TYPE( POP_Native )         :: irfFields
   CHARACTER(200)             :: thisIRFFile, meshFile
   CHARACTER(400)             :: crsFile
   CHARACTER(17)              :: walltime
   CHARACTER(5)               :: fileIDChar
   INTEGER                    :: nEntries, irf_id, row, col, m, diff_id
   INTEGER                    :: i, j, k, this_i, this_j, this_k
   INTEGER                    :: true_i, true_j, kdiff
   INTEGER                    :: fUnit, fileID
   INTEGER, ALLOCATABLE       :: nval(:), columns(:,:)
   REAL(prec), ALLOCATABLE    :: opdata(:,:)
   INTEGER, ALLOCATABLE       :: nval_diff(:), columns_diff(:,:)
   REAL(prec), ALLOCATABLE    :: opdata_diff(:,:)
   REAL(prec)                 :: t0, t1    


      CALL params % Build( )

      ! Use the first IRF file (assumed to list netcdf files) to obtain the
      ! meshfile
      PRINT*, TRIM(params % IRFListFile)
      OPEN( UNIT=NewUnit(fUnit),&
            FILE=TRIM(params % IRFListFile), &
            FORM='FORMATTED',&
            ACCESS='SEQUENTIAL',&
            ACTION='READ',&
            STATUS='OLD' )
 
      READ( fUnit, '(A200)' ) meshFile
      CLOSE( fUnit )
      PRINT*,meshFile
      ! Load in the mesh from the netcdf file specified above
      CALL mesh % Load( TRIM(meshfile) )
      mesh % meshType = params % meshType ! And set a flag for the type of mesh

      ! Here, we build a stencil for the finite difference scheme
      ! specified above, without overlaps
      CALL advStencil % Build( stencilFlag = params % StencilType, &
                               flavor      = Normal )

      ! Allocate space for the transport operator
      ! The number of possible non-zero entries is the number of degrees of
      ! freedom multiplied by the stencil-width
      nEntries = ( mesh % nDOF )*( advStencil % nPoints )
      CALL transportOperator % Build( mesh % nDOF, mesh % nDOF, nEntries )
      CALL diffusionOperator % Build( mesh % nDOF, mesh % nDOF, mesh % nDOF*3 ) 

      ! Allocate space for a hash-table storage
      ALLOCATE( nval(1:mesh % nDOF), columns(1:mesh % nDOF, 1:advStencil % nPoints) )
      ALLOCATE( opdata(1:mesh % nDOF, 1:advStencil % nPoints) )
      ! Allocate space for a hash-table storage of vertical diffusion operator
      ALLOCATE( nval_diff(1:mesh % nDOF), columns_diff(1:mesh % nDOF, 1:3) )
      ALLOCATE( opdata_diff(1:mesh % nDOF, 1:3) )


      ! Read the adjacency graph from file
      CALL graph % ReadGraphBinFile( TRIM(params % graphFile) )
      
      ! Allocate space for the IRF fields. The number of IRF fields is
      ! equivalent to the number of colors in the graph
      !
      ! When building this data structure for the irf-fields, nColors
      ! corresponds to the number of irf fields-we add an additionation
      ! "POP_Native" attribute for the vertical diffusivity
      CALL irfFields % Build( mesh, graph % nColors+1 )
      diff_id = graph % nColors + 1
      
      ! /////////////////////// Operator Diagnosis //////////////////////////// !

      OPEN( UNIT=NewUnit(fUnit),&
            FILE=TRIM(params % IRFListFile), &
            FORM='FORMATTED',&
            ACCESS='SEQUENTIAL',&
            ACTION='READ',&
            STATUS='OLD' )
 
      DO fileID = 1, params % nIRFFiles
 
         READ( fUnit, '(A200)' ) thisIRFFile

         IF( fileID >= params % IRFStart )THEN

         ! Initialize the NetCDF file
         PRINT*, 'Initialize for reading...'//TRIM( thisIRFFile )
         CALL irfFields % InitializeForNetCDFRead( modelType=ImpulseResponseField, &
                                                   filename=TRIM(thisIRFFile), &
                                                   initOn=.TRUE. )
  
         ! Read in all of the impulse response fields
         PRINT*, 'Loading Impulse Response Functions.'
         CALL irfFields % LoadTracerFromNetCDF( mesh ) 
         ! Close the netcdf file
         CALL irfFields % FinalizeNetCDF( )

         PRINT*, 'Diagnosing Transport Operators.'
         CALL CPU_TIME( t0 )
         nval    = 0
         columns = 0
         opdata  = 0.0_prec
         DO irf_id = 1, graph % nColors
            DO col = 1, mesh % nDOF

               IF( graph % color(col) == irf_id )THEN

                  i = mesh % dofToIJK(1,col)
                  j = mesh % dofToIJK(2,col)
                  k = mesh % dofToIJK(3,col)

                  DO m = 1, advstencil % nPoints

                     this_i = i + advstencil % relativeNeighbors(1,m)
                     this_j = j + advstencil % relativeNeighbors(2,m)
                     this_k = k + advstencil % relativeNeighbors(3,m)
                 
                     IF( this_k <= mesh % nZ .AND. this_k > 0 )THEN

                        
                        CALL GetTrueIJ( params % MeshType, &
                                        this_i, this_j, &
                                        mesh % nX, mesh % nY, &
                                        true_i, true_j )
                           IF( true_i > mesh % nX .OR. true_j > mesh % nY )THEN
                              PRINT*, 'Code attempting DOF data access out of bounds! STOPPING!'
                              STOP
                           ENDIF

                        row = mesh % ijkToDOF(true_i,true_j,this_k)

                        IF( row /= 0 )THEN
                           IF( row > mesh % nDOF )THEN
                              PRINT*, 'Code attempting ROW data access out of bounds! STOPPING!'
                              STOP
                           ENDIF

                           nval(row) = nval(row) + 1

                           IF( nVal(row) > advStencil % nPoints )THEN
                              PRINT*, 'Code attempting COLUMN data access out of bounds! STOPPING!'
                              STOP
                           ENDIF

                           opdata(row,nval(row))  = irfFields % tracer(true_i,true_j,this_k,irf_id)
                           columns(row,nval(row)) = col
                        ENDIF

                     ENDIF

                  ENDDO

               ENDIF

            ENDDO
         ENDDO

         PRINT*, 'Diagnosing Vertical Diffusion Operator.'
         ! Diffusion operator diagnosis
         nval_diff    = 0
         columns_diff = 0
         opdata_diff  = 0.0_prec
         DO row = 1, mesh % nDOF

            i = mesh % dofToIJK(1,row)
            j = mesh % dofToIJK(2,row)
            k = mesh % dofToIJK(3,row)

            IF( k <= mesh % KMT(i,j) )THEN

            IF( k == 1 )THEN ! At the surface, the "no-diffusive-flux" condition is used
                ! cell k coefficient/diagonal
                col = row
                nval_diff(row)                   = nval_diff(row) + 1
                opdata_diff(row,nval_diff(row))  = -irfFields % tracer(i,j,k,diff_id)/mesh % dz(k)/mesh % dzw(k)
                columns_diff(row,nval_diff(row)) =  col

                ! cell k+1
                col = mesh % ijkToDOF(i,j,k+1)
                nval_diff(row)                   = nval_diff(row) + 1
                opdata_diff(row,nval_diff(row))  = irfFields % tracer(i,j,k,diff_id)/mesh % dz(k)/mesh % dzw(k)
                columns_diff(row,nval_diff(row)) =  col

            ELSEIF( k == mesh % KMT(i,j) )THEN !At the bottom-- no diffusive flux

                ! cell k-1
                col = mesh % ijkToDOF(i,j,k-1)
                nval_diff(row)                   = nval_diff(row) + 1
                opdata_diff(row,nval_diff(row))  = irfFields % tracer(i,j,k-1,diff_id)/mesh % dz(k)/mesh % dzw(k-1)
                columns_diff(row,nval_diff(row)) = col

                ! cell k coefficient/diagonal
                col = row
                nval_diff(row)                   = nval_diff(row) + 1
                opdata_diff(row,nval_diff(row))  = -irfFields % tracer(i,j,k-1,diff_id)/mesh % dz(k)/mesh % dzw(k-1)
                columns_diff(row,nval_diff(row)) =  col

            ELSE

                ! cell k-1
                col = mesh % ijkToDOF(i,j,k-1)
                nval_diff(row)                   = nval_diff(row) + 1
                opdata_diff(row,nval_diff(row))  = irfFields % tracer(i,j,k-1,diff_id)/mesh % dz(k)/mesh % dzw(k-1)
                columns_diff(row,nval_diff(row)) = col

                ! cell k coefficient/diagonal
                col = row
                nval_diff(row)                   = nval_diff(row) + 1
                opdata_diff(row,nval_diff(row))  = -irfFields % tracer(i,j,k-1,diff_id)/mesh % dz(k)/mesh % dzw(k-1)-&
                                                    irfFields % tracer(i,j,k,diff_id)/mesh % dz(k)/mesh % dzw(k)
                columns_diff(row,nval_diff(row)) =  col

                ! cell k+1
                col = mesh % ijkToDOF(i,j,k+1)
                nval_diff(row)                   = nval_diff(row) + 1
                opdata_diff(row,nval_diff(row))  = irfFields % tracer(i,j,k,diff_id)/mesh % dz(k)/mesh % dzw(k)
                columns_diff(row,nval_diff(row)) = col

            ENDIF
            ENDIF                 
         ENDDO
         CALL CPU_TIME( t1 )
         PRINT*, 'DONE! '
         WRITE( walltime, '(F17.4)' ) t1-t0
         PRINT*, 'Time to diagnose operator : '//TRIM(walltime)//' sec'

         PRINT*, 'Converting hash-table...'

         ! Convert from the "hash-table" to the CRS Format
         kdiff = 1
         k     = 1
         DO row = 1, mesh % nDOF
            IF( nval(row) > 0 )THEN

               transportOperator % rowBounds(1,row) = k
               DO m = 1, nval(row)
                  transportOperator % A(k)   = opdata(row,m)
                  transportOperator % col(k) = columns(row,m)
                  k = k+1
               ENDDO
               transportOperator % rowBounds(2,row) = k-1

            ENDIF
            IF( nval_diff(row) > 0 )THEN

               diffusionOperator % rowBounds(1,row) = kdiff
               DO m = 1, nval_diff(row)
                  diffusionOperator % A(kdiff)   = opdata_diff(row,m)
                  diffusionOperator % col(kdiff) = columns_diff(row,m)
                  kdiff = kdiff+1
               ENDDO
               diffusionOperator % rowBounds(2,row) = kdiff-1
            ENDIF
         ENDDO
      
         WRITE(fileIDChar, '(I5.5)' ) fileID
         crsFile=TRIM(params % feotsOperatorDirectory)//TRIM(params % operatorBaseName)//'_advect.'//fileIDChar
         PRINT*,'Writing CRS Matrix files : '//TRIM(crsFile)
         CALL transportOperator % WriteSparseConnectivity( TRIM(crsFile) )
         CALL transportOperator % WriteMatrixData( TRIM(crsFile) )

         crsFile=TRIM(params % feotsOperatorDirectory)//TRIM(params % operatorBaseName)//'_vdiffu.'//fileIDChar
         PRINT*,'Writing CRS Matrix files : '//TRIM(crsFile)
         CALL diffusionOperator % WriteSparseConnectivity( TRIM(crsFile) )
         CALL diffusionOperator % WriteMatrixData( TRIM(crsFile) )


         CALL transportOperator % Reset( )
         CALL diffusionOperator % Reset( )

         ENDIF

      ENDDO ! Loop over the IRF Files
      CLOSE( fUnit )

      ! Clear memory
      CALL advStencil % Trash( )
      CALL graph % Trash( )
      CALL mesh % Trash( )
      CALL irfFields % Trash( )
      CALL transportOperator % Trash( )
      CALL diffusionOperator % Trash( )
      DEALLOCATE( nval, columns )
      DEALLOCATE( opdata )

  END SUBROUTINE OperatorDiagnosis

  SUBROUTINE RegionalExtraction( )

    IMPLICIT NONE

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

      CALL params % Build( )

      CALL globalMesh % Load( TRIM( params % meshFile ) )

      CALL modelstencil % Build( stencilFlag = params % stencilType, &
                                 flavor      = LateralPlusCorners )

      CALL region % Build( globalMesh, regionalMesh, modelstencil, params % meshType, &
                           params % south, params % north, &
                           params % east, params % west, params % maskfile ) 

      CALL regionalMesh % WriteNetCDF( TRIM(params % regionalMeshFile) )
      ! Write the regional data structure to a pickup file for later use

      IF( TRIM(params % maskfile) == '' )THEN
         CALL region % WritePickup( TRIM(params % regionalOperatorDirectory)//'mappings', maskProvided=.FALSE. )
      ELSE
         CALL region % WritePickup( TRIM(params % regionalOperatorDirectory)//'mappings', maskProvided=.TRUE. )
      ENDIF

      IF( params % ExtractRegionalOperators )THEN

         CALL advstencil % Build( stencilFlag = params % stencilType, &
                                  flavor      = Normal )

         nGentries = ( globalMesh % nDOF )*( advStencil % nPoints )
         CALL transportOp % Build( globalmesh % nDOF, globalmesh % nDOF, nGEntries )
         CALL diffusionOp % Build( globalmesh % nDOF, globalmesh % nDOF, globalmesh % nDOF*3 ) 
         
         nRentries = ( region % nCells )*( advStencil % nPoints )
         CALL regionalTransportOp % Build( region % nCells, region % nCells, nREntries )
         CALL regionalDiffusionOp % Build( region % nCells, region % nCells, region % nCells*3 ) 

         PRINT*, ' Extracting regional operators.'
         DO fileID = 1, params % nIRFFiles

            IF( fileID >= params % IRFStart )THEN

               WRITE(fileIDChar, '(I5.5)' ) fileID
               crsFile=TRIM(params % feotsOperatorDirectory)//TRIM(params % operatorBaseName)//'_advect.'//fileIDChar
               PRINT*,'Reading CRS Matrix files : '//TRIM(crsFile)
               CALL transportOp % ReadSparseConnectivity( TRIM(crsFile) )
               CALL transportOp % ReadMatrixData( TRIM(crsFile) )

               crsFile=TRIM(params % feotsOperatorDirectory)//TRIM(params % operatorBaseName)//'_vdiffu.'//fileIDChar
               PRINT*,'Reading CRS Matrix files : '//TRIM(crsFile)
               CALL diffusionOp % ReadSparseConnectivity( TRIM(crsFile) )
               CALL diffusionOp % ReadMatrixData( TRIM(crsFile) )

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

               crsFile=TRIM(params % regionalOperatorDirectory)//TRIM(params % operatorBaseName)//'_advect.'//fileIDChar
               PRINT*,'Writing CRS Matrix files : '//TRIM(crsFile)
               CALL regionalTransportOp % WriteSparseConnectivity( TRIM(crsFile) )
               CALL regionalTransportOp % WriteMatrixData( TRIM(crsFile) )

               crsFile=TRIM(params % regionalOperatorDirectory)//TRIM(params % operatorBaseName)//'_vdiffu.'//fileIDChar
               PRINT*,'Writing CRS Matrix files : '//TRIM(crsFile)
               CALL regionalDiffusionOp % WriteSparseConnectivity( TRIM(crsFile) )
               CALL regionalDiffusionOp % WriteMatrixData( TRIM(crsFile) )
            
            ENDIF
         ENDDO
   
         CALL transportOp % Trash( )
         CALL diffusionOp % Trash( )
         CALL regionaltransportOp % Trash( )
         CALL regionalDiffusionOp % Trash( )
        
      ENDIF   
      
      ! Clean up memory !
     ! DEALLOCATE( maskfield )
      CALL globalMesh % Trash( )
      CALL regionalMesh % Trash( )
      CALL modelstencil % Trash( )
      CALL region % Trash( )    

  END SUBROUTINE RegionalExtraction

  SUBROUTINE FEOTSInitialize()
    
    IMPLICIT NONE

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
      !CALL feots % Build( myRank == 0, nProcs == 1 )
      CALL feots % Build( myRank, nProcs )

      CALL InitialConditions( feots )

      !  //////////////////////////////////////////// File I/O  //////////////////////////////////////////////////////// !
      CALL feots % nativeSol % InitializeForNetCDFWrite( feots % params % TracerModel, &
                                                         feots % mesh, &
                                                         TRIM(feots % params % outputDirectory)//'Tracer.init.nc', &
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

 SUBROUTINE FEOTSIntegrate()

   IMPLICIT NONE

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

      CALL feots % Build( myRank, nProcs )

      recordID = 1

         IF( myRank == 0 )THEN 
            fileID   = feots % params % iterInit
            ! /////////////////////////// Load in the initial conditions //////////////////// !
            WRITE( ncfileTag, '(I10.10)' ) fileID
            ! For now, the record ID is 1. In the future, this will need to be
            ! calculated as a function of the initial iterate, the dump frequency,
            ! and the number of records per netcdf file
            
           ! Tracer.init.nc is read for the mask and source terms
            CALL feots % nativeSol % InitializeForNetCDFRead( feots % params % TracerModel,'Tracer.init.nc', .TRUE. )
            CALL feots % nativeSol % ReadSourceEtcNetCDF( feots % mesh )
            CALL feots % nativeSol % ReadNetCDFRecord( feots % mesh, recordID )
            CALL feots % nativeSol % FinalizeNetCDF( )


            ! ************
            ! This section of code loads in the temperature, salinity, potential
            ! density, and ssh fields if the water mass tagging is turned on.
            ! In future implementations with the volume correction, this will need
            ! to be turned on the volume corrections are enabled. 
            !
            IF( feots % params % WaterMassTagging ) THEN

               WRITE( fileIDChar, '(I5.5)' ) feots % params % IRFStart
               IF( feots % params % Regional )THEN
                  CALL feots % nativeSol % LoadOceanState( feots % mesh, &
                                                          TRIM(feots % params % regionalOperatorDirectory)//'Ocean.'//fileIDChar//'.nc')
               ELSE
                  OPEN( UNIT=NewUnit(fUnit),&
                        FILE=TRIM(feots % params % IRFListFile), &
                        FORM='FORMATTED',&
                        ACCESS='SEQUENTIAL',&
                        ACTION='READ',&
                        STATUS='OLD' )
            
                  DO fileID = 1, feots % params % nIRFFiles
            
                     READ( fUnit, '(A200)' ) thisIRFFile
            
                     IF( fileID == feots % params % IRFStart )THEN
                        CALL feots % nativeSol % LoadOceanState( feots % mesh,TRIM(thisIRFFile) )
                     ENDIF
                  ENDDO
   
                  CLOSE(fUnit)
               ENDIF
            ENDIF

         ENDIF ! myRank == 0
         !***********

         IF( feots % params % iterInit == 0 )THEN
            tn = 0.0_prec
         ELSE

            IF( myRank == 0 )THEN 
               ! This pickup file is read for the correct "initial condition"
               CALL feots % nativeSol % InitializeForNetCDFRead( feots % params % TracerModel,'Tracer.'//ncFileTag//'.nc', .FALSE. )
               CALL feots % nativeSol % ReadNetCDFRecord( feots % mesh, recordID )
               CALL feots % nativeSol % FinalizeNetCDF( )
            ENDIF
               tn = REAL(feots % params % iterInit,prec)*feots % params % dt

         ENDIF
   
         ! /////////////////////////////////////////////////////////////////////////////// !
         
         ! Transfer the data from the native storage to the FEOTS storage
         IF( myRank == 0 )THEN
            CALL feots % MapAllToDOF( )
         ENDIF
#ifdef HAVE_MPI
         CALL MPI_BARRIER( MPI_COMM_WORLD, mpiErr )
         CALL feots % ScatterSolution( myRank, nProcs )
         CALL feots % ScatterSource( myRank, nProcs )
         CALL feots % ScatterMask( myRank, nProcs )
         IF( myRank /= 0 )THEN
            CALL feots % MapTracerToDOF( )
         ENDIF
#endif 
         ! //// Forward Mode //// !
         PRINT*, '  Starting ForwardStep'
   
         DO iter = feots % params % iterInit, feots % params % iterInit + feots % params % nTimeSteps -1, feots % params % nStepsPerDump

            !$OMP PARALLEL
            IF( myRank /= 0 .OR. nProcs == 0)THEN
               CALL feots % ForwardStep( tn, feots % params % nStepsPerDump, myRank, nProcs )
            ENDIF
            !$OMP END PARALLEL 
 

            IF(myRank /= 0  .OR. nProcs == 0)THEN
               CALL feots % MapTracerFromDOF( )
            ENDIF

#ifdef HAVE_MPI
            CALL feots % GatherSolution( myRank, nProcs )
#endif
            IF( myRank == 0 )THEN

               WRITE( ncfileTag, '(I10.10)' ) iter + feots % params % nStepsPerDump
               CALL feots % nativeSol % InitializeForNetCDFWrite( feots % params % TracerModel, &
                                                                  feots % mesh, &
                                                                 'Tracer.'//ncFileTag//'.nc', &
                                                                  .FALSE. )
               CALL feots % nativeSol % WriteNetCDFRecord( feots % mesh, recordID )
               CALL feots % nativeSol % FinalizeNetCDF( )
            ENDIF

#ifdef HAVE_MPI
            CALL MPI_BARRIER( MPI_COMM_WORLD, mpiErr )
#endif

         ENDDO

      CALL feots % Trash( )
#ifdef HAVE_MPI
      CALL MPI_FINALIZE( mpiErr )
#endif

 END SUBROUTINE FEOTSIntegrate

 SUBROUTINE FEOTSEquilibrate()

   IMPLICIT NONE

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

      CALL feots % Build( myRank, nProcs )

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
#ifdef HAVE_MPI
      CALL MPI_BARRIER( MPI_COMM_WORLD, mpiErr )
      CALL feots % ScatterSolution( myRank, nProcs )
      CALL feots % ScatterSource( myRank, nProcs )
      CALL feots % ScatterMask( myRank, nProcs )
#endif 

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
