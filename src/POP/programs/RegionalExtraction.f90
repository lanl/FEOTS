PROGRAM RegionalExtraction

! src/matrices/
USE CRSMatrix_Class
! src/POP/
USE POP_Mesh_Class
USE POP_Stencil_Class
USE POP_Regional_Class
USE POP_Params_Class


IMPLICIT NONE

   TYPE( CRSMatrix )    :: transportOp, diffusionOP
   TYPE( CRSMatrix )    :: regionalTransportOp, regionalDiffusionOP
   TYPE( POP_Params )   :: params
   TYPE( POP_Mesh )     :: globalMesh 
   TYPE( POP_Mesh )     :: regionalMesh 
   TYPE( Stencil )      :: modelstencil, advStencil
   TYPE( POP_Regional ) :: region
   
   INTEGER        :: fileID, nGentries, nRentries
   INTEGER, ALLOCATABLE :: maskfield(:,:)
   CHARACTER(5)   :: fileIDChar
   CHARACTER(400) :: crsFile

      CALL params % Build( )

      CALL globalMesh % Load( TRIM( params % meshFile ) )

      ALLOCATE( maskfield(1:globalmesh % nX, 1:globalmesh % nY) )
      CALL modelstencil % Build( stencilFlag = params % stencilType, &
                                 flavor      = LateralPlusCorners )

      IF( TRIM(params % maskfile) == '' )THEN
         CALL region % Build( globalMesh, modelstencil, params % meshType, &
                              params % south, params % north, &
                              params % east, params % west ) 
         CALL region % GenerateRegionalMesh( globalMesh, regionalMesh )
      ELSE

         PRINT*, ' Using Mask file ', TRIM(params % maskFile)
         CALL region % LoadMaskField( globalmesh, maskfield, params % maskfile )

         CALL region % Build( globalMesh, modelstencil, params % meshType, &
                              params % south, params % north, &
                              params % east, params % west, maskfield ) 
         CALL region % GenerateRegionalMesh( globalMesh, regionalMesh, maskfield )
      ENDIF

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
      DEALLOCATE( maskfield )
      CALL globalMesh % Trash( )
      CALL regionalMesh % Trash( )
      CALL modelstencil % Trash( )
      CALL region % Trash( )    

END PROGRAM RegionalExtraction
