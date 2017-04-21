PROGRAM OperatorDiagnosis

! src/common/
USE ModelPrecision
USE CommonRoutines
! src/matrices/
USE CRSMatrix_Class
! src/POP/
USE POP_Params_Class
USE POP_Mesh_Class
USE POP_Stencil_Class
USE POP_AdjacencyGraph_Class
USE POP_Native_Class
USE POP_GridTypeMappings

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
                                                   filename=TRIM(thisIRFFile) )
  
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
         DO col = 1, mesh % nDOF

            i = mesh % dofToIJK(1,col)
            j = mesh % dofToIJK(2,col)
            k = mesh % dofToIJK(3,col)

            IF( k <= mesh % KMT(i,j) )THEN

            IF( k == 1 )THEN ! At the surface, the "no-diffusive-flux" condition is used
                ! cell k coefficient/diagonal
                row = col
                nval_diff(row)                   = nval_diff(row) + 1
                opdata_diff(row,nval_diff(row))  = -irfFields % tracer(i,j,k,diff_id)/mesh % dz(k)/mesh % dzw(k)
                columns_diff(row,nval_diff(row)) =  col

                ! cell k+1
                row = mesh % ijkToDOF(i,j,k+1)
                nval_diff(row)                   = nval_diff(row) + 1
                opdata_diff(row,nval_diff(row))  = irfFields % tracer(i,j,k,diff_id)/mesh % dz(k)/mesh % dzw(k)
                columns_diff(row,nval_diff(row)) =  col

            ELSEIF( k == mesh % KMT(i,j) )THEN !At the bottom-- no diffusive flux

                ! cell k-1
                row = mesh % ijkToDOF(i,j,k-1)
                nval_diff(row)                   = nval_diff(row) + 1
                opdata_diff(row,nval_diff(row))  = irfFields % tracer(i,j,k-1,diff_id)/mesh % dz(k)/mesh % dzw(k-1)
                columns_diff(row,nval_diff(row)) =  col

                ! cell k coefficient/diagonal
                row = col
                nval_diff(row)                   = nval_diff(row) + 1
                opdata_diff(row,nval_diff(row))  = -irfFields % tracer(i,j,k-1,diff_id)/mesh % dz(k)/mesh % dzw(k-1)
                columns_diff(row,nval_diff(row)) =  col

            ELSE

                ! cell k-1
                row = mesh % ijkToDOF(i,j,k-1)
                nval_diff(row)                   = nval_diff(row) + 1
                opdata_diff(row,nval_diff(row))  = irfFields % tracer(i,j,k,diff_id)/mesh % dz(k)/mesh % dzw(k-1)
                columns_diff(row,nval_diff(row)) =  col

                ! cell k coefficient/diagonal
                row = col
                nval_diff(row)                   = nval_diff(row) + 1
                opdata_diff(row,nval_diff(row))  = -irfFields % tracer(i,j,k-1,diff_id)/mesh % dz(k)/mesh % dzw(k-1)-&
                                                    irfFields % tracer(i,j,k,diff_id)/mesh % dz(k)/mesh % dzw(k)
                columns_diff(row,nval_diff(row)) =  col

                ! cell k+1
                row = mesh % ijkToDOF(i,j,k+1)
                nval_diff(row)                   = nval_diff(row) + 1
                opdata_diff(row,nval_diff(row))  = irfFields % tracer(i,j,k,diff_id)/mesh % dz(k)/mesh % dzw(k)
                columns_diff(row,nval_diff(row)) =  col

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

END PROGRAM OperatorDiagnosis

