! ExtractOceanState.f90
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

PROGRAM ExtractOceanState

! src/POP/
USE POP_FEOTS_Class
USE POP_Mesh_Class
USE POP_Native_Class


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

END PROGRAM ExtractOceanState
