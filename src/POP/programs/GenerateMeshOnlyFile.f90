! GenerateMeshOnlyFile.f90
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


PROGRAM GenerateMeshOnlyFile

! src/common/
USE CommonRoutines
! src/POP/
USE POP_Params_Class
USE POP_Mesh_Class

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

END PROGRAM GenerateMeshOnlyFile

