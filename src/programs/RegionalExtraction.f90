! RegionalExtraction.f90
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

END PROGRAM RegionalExtraction
