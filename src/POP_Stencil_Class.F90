! POP_Stencil_Class.f90
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

 
MODULE POP_Stencil_Class


! src/common/
USE ConstantsDictionary

IMPLICIT NONE

   TYPE Stencil
      INTEGER              :: nPoints ! Number of points in the stencil
      INTEGER, ALLOCATABLE :: relativeNeighbors(:,:) ! indices of the neighbors relative to the center

      CONTAINS 
  
      PROCEDURE :: Build => Build_Stencil
      PROCEDURE :: Trash => Trash_Stencil

   END TYPE Stencil


CONTAINS

 SUBROUTINE Build_Stencil( myStencil, stencilFlag, flavor )
 ! This routine handles the construction of the Stencil class based on pre-defined 
 ! stencils.
 ! The stencil flag refers to a particular advection or diffusion stencil
 !   Currently, stencilFlag == LaxWendroff is only supported.
 ! The stencil can come in many flavors.
 !   Flavor == normal  refers to the finite difference stencil
 !   Flavor == overlap refers to the the stencil obtained by applying the finite difference stencil
 !                     to each neighbor in the normal stencil.
 !   Flavor == lateral refers to the finite difference stencil without the vertical neighbors
   IMPLICIT NONE
   CLASS( Stencil ), INTENT(inout) :: myStencil
   INTEGER, INTENT(in)             :: stencilFlag
   INTEGER, INTENT(in)             :: flavor
   ! Local
   INTEGER :: i, j, k, iStencil

      IF( stencilFlag == LaxWendroff .AND. flavor == Overlap )THEN
         
         PRINT*, 'S/R : Build_Stencil : Constructing Lax-Wendroff with overlap flavor.'
         myStencil % nPoints = 45
         ALLOCATE( myStencil % relativeNeighbors(1:3,1:45) )

         iStencil = 1

         k = -2 ! Vertical Level
         myStencil % relativeNeighbors(1:3,iStencil) = (/ 0,0,k /)

         k = -1 ! Vertical Level
         DO j = -1, 1
           DO i = -1, 1

             iStencil = iStencil + 1
             myStencil % relativeNeighbors(1:3,iStencil) = (/ i,j,k /)
            
           ENDDO
         ENDDO

         k = 0 ! Vertical Level
         DO j = -2, 2
           DO i = -2, 2

             iStencil = iStencil + 1
             myStencil % relativeNeighbors(1:3,iStencil) = (/ i,j,k /)
            
           ENDDO
         ENDDO

         k = 1 ! Vertical Level
         DO j = -1, 1
           DO i = -1, 1

             iStencil = iStencil + 1
             myStencil % relativeNeighbors(1:3,iStencil) = (/ i,j,k /)
            
           ENDDO
         ENDDO

         k = 2 ! Vertical Level
         iStencil = iStencil + 1
         myStencil % relativeNeighbors(1:3,iStencil) = (/ 0,0,k /)

      ELSEIF( stencilFlag == LaxWendroff .AND. flavor == Normal )THEN
         
         PRINT*, 'S/R : Build_Stencil : Constructing Lax-Wendroff stencil with Normal flavor.'
         myStencil % nPoints = 11
         ALLOCATE( myStencil % relativeNeighbors(1:3,1:11) )

         iStencil = 1

         k = -1 ! Vertical Level
         myStencil % relativeNeighbors(1:3,iStencil) = (/ 0,0,k /)

         k = 0 ! Vertical Level
         DO j = -1, 1
           DO i = -1, 1

             iStencil = iStencil + 1
             myStencil % relativeNeighbors(1:3,iStencil) = (/ i,j,k /)
            
           ENDDO
         ENDDO

         k = 1 ! Vertical Level
         iStencil = iStencil + 1
         myStencil % relativeNeighbors(1:3,iStencil) = (/ 0,0,k /)

      ELSEIF( stencilFlag == LaxWendroff27 .AND. flavor == Overlap )THEN
         
         PRINT*, 'S/R : Build_Stencil : Constructing Lax-Wendroff 27-point stencil with Overlap flavor.'
         myStencil % nPoints = 125
         ALLOCATE( myStencil % relativeNeighbors(1:3,1:125) )

         iStencil = 0

         DO k = -2, 2 ! Vertical Level
           DO j = -2, 2
             DO i = -2, 2

               iStencil = iStencil + 1
               myStencil % relativeNeighbors(1:3,iStencil) = (/ i,j,k /)
              
             ENDDO
           ENDDO
         ENDDO

      ELSEIF( stencilFlag == LaxWendroff27 .AND. flavor == Normal )THEN
         
         PRINT*, 'S/R : Build_Stencil : Constructing Lax-Wendroff 27-point stencil with Normal flavor.'
         myStencil % nPoints = 27
         ALLOCATE( myStencil % relativeNeighbors(1:3,1:27) )

         iStencil = 0

         DO k = -1, 1 ! Vertical Level
           DO j = -1, 1
             DO i = -1, 1

               iStencil = iStencil + 1
               myStencil % relativeNeighbors(1:3,iStencil) = (/ i,j,k /)
              
             ENDDO
           ENDDO
         ENDDO

      ELSEIF( stencilFlag == Upwind3 .AND. flavor == Overlap )THEN
         
         PRINT*, 'S/R : Build_Stencil : Constructing 3rd-order upwind stencil with Overlap flavor.'
         myStencil % nPoints = 73
         ALLOCATE( myStencil % relativeNeighbors(1:3,1:73) )

         iStencil = 0

         DO k = -4, -3 ! Vertical Level
           iStencil = iStencil + 1
           myStencil % relativeNeighbors(1:3,iStencil) = (/ 0,0,k /)
         ENDDO

         DO k = -2, -1 ! Vertical Level
           DO j = -2, 2
             iStencil = iStencil + 1
             myStencil % relativeNeighbors(1:3,iStencil) = (/ 0,j,k /)
           ENDDO
           DO i = -2, -1
             iStencil = iStencil + 1
             myStencil % relativeNeighbors(1:3,iStencil) = (/ i,0,k /)
           ENDDO
           DO i = 1, 2
             iStencil = iStencil + 1
             myStencil % relativeNeighbors(1:3,iStencil) = (/ i,0,k /)
           ENDDO
         ENDDO

         k = 0 ! Vertical Level
         DO j = -4, -3
           iStencil = iStencil + 1
           myStencil % relativeNeighbors(1:3,iStencil) = (/ 0,j,k /)
         ENDDO
         DO i = -4, -3
           iStencil = iStencil + 1
           myStencil % relativeNeighbors(1:3,iStencil) = (/ i,0,k /)
         ENDDO
         DO j = -2, 2
           DO i = -2, 2
             iStencil = iStencil + 1
             myStencil % relativeNeighbors(1:3,iStencil) = (/ i,j,k /)
           ENDDO
         ENDDO
         DO j = 3, 4
           iStencil = iStencil + 1
           myStencil % relativeNeighbors(1:3,iStencil) = (/ 0,j,k /)
         ENDDO
         DO i = 3, 4
           iStencil = iStencil + 1
           myStencil % relativeNeighbors(1:3,iStencil) = (/ i,0,k /)
         ENDDO

         DO k = 1, 2 ! Vertical Level
           DO j = -2, 2
             iStencil = iStencil + 1
             myStencil % relativeNeighbors(1:3,iStencil) = (/ 0,j,k /)
           ENDDO
           DO i = -2, -1
             iStencil = iStencil + 1
             myStencil % relativeNeighbors(1:3,iStencil) = (/ i,0,k /)
           ENDDO
           DO i = 1, 2
             iStencil = iStencil + 1
             myStencil % relativeNeighbors(1:3,iStencil) = (/ i,0,k /)
           ENDDO
         ENDDO

         DO k = 3, 4 ! Vertical Level
           iStencil = iStencil + 1
           myStencil % relativeNeighbors(1:3,iStencil) = (/ 0,0,k /)
         ENDDO

      ELSEIF( stencilFlag == Upwind3 .AND. flavor == Normal )THEN
         
         PRINT*, 'S/R : Build_Stencil : Constructing 3rd-order upwind stencil with Normal flavor.'
         myStencil % nPoints = 13
         ALLOCATE( myStencil % relativeNeighbors(1:3,1:13) )

         iStencil = 0

         DO k = -2, -1 ! Vertical Level
           iStencil = iStencil + 1
           myStencil % relativeNeighbors(1:3,iStencil) = (/ 0,0,k /)
         ENDDO

         k = 0 ! Vertical Level
         DO j = -2, 2
           iStencil = iStencil + 1
           myStencil % relativeNeighbors(1:3,iStencil) = (/ 0,j,k /)
         ENDDO
         DO i = -2, -1
           iStencil = iStencil + 1
           myStencil % relativeNeighbors(1:3,iStencil) = (/ i,0,k /)
         ENDDO
         DO i = 1, 2
           iStencil = iStencil + 1
           myStencil % relativeNeighbors(1:3,iStencil) = (/ i,0,k /)
         ENDDO

         DO k = 1, 2 ! Vertical Level
           iStencil = iStencil + 1
           myStencil % relativeNeighbors(1:3,iStencil) = (/ 0,0,k /)
         ENDDO

      ELSEIF( stencilFlag == LaxWendroff .AND. flavor == Lateral )THEN
         
         PRINT*, 'S/R : Build_Stencil : Constructing Lax-Wendroff stencil with Lateral flavor.'
         myStencil % nPoints = 9
         ALLOCATE( myStencil % relativeNeighbors(1:3,1:9) )

         iStencil = 0
         k = 0 ! Vertical Level
         DO j = -1, 1
           DO i = -1, 1

             iStencil = iStencil + 1
             myStencil % relativeNeighbors(1:3,iStencil) = (/ i,j,k /)
            
           ENDDO
         ENDDO

      ELSEIF( stencilFlag == LaxWendroff27 .AND. flavor == Lateral )THEN
         
         PRINT*, 'S/R : Build_Stencil : Constructing Lax-Wendroff 27 point stencil with Lateral flavor.'
         myStencil % nPoints = 9
         ALLOCATE( myStencil % relativeNeighbors(1:3,1:9) )

         iStencil = 0
         k = 0 ! Vertical Level
         DO j = -1, 1
           DO i = -1, 1

             iStencil = iStencil + 1
             myStencil % relativeNeighbors(1:3,iStencil) = (/ i,j,k /)
            
           ENDDO
         ENDDO
      ELSEIF( stencilFlag == LaxWendroff .AND. flavor == LateralPlusCorners )THEN
         
         PRINT*, 'S/R : Build_Stencil : Constructing Lax-Wendroff stencil with LateralPlusCorners flavor.'
         myStencil % nPoints = 9
         ALLOCATE( myStencil % relativeNeighbors(1:3,1:9) )

         iStencil = 0
         k = 0 ! Vertical Level
         DO j = -1, 1
           DO i = -1, 1

             iStencil = iStencil + 1
             myStencil % relativeNeighbors(1:3,iStencil) = (/ i,j,k /)
            
           ENDDO
         ENDDO

      ELSEIF( stencilFlag == LaxWendroff27 .AND. flavor == LateralPlusCorners )THEN
         
         PRINT*, 'S/R : Build_Stencil : Constructing Lax-Wendroff 27-point stencil with LateralPlusCorners flavor.'
         myStencil % nPoints = 9
         ALLOCATE( myStencil % relativeNeighbors(1:3,1:9) )

         iStencil = 0
         k = 0 ! Vertical Level
         DO j = -1, 1
           DO i = -1, 1

             iStencil = iStencil + 1
             myStencil % relativeNeighbors(1:3,iStencil) = (/ i,j,k /)
            
           ENDDO
         ENDDO

      ELSE
         PRINT*, 'S/R : Build_Stencil : Invalid Stencil Flag! Stopping'
         STOP
      ENDIF

 END SUBROUTINE Build_Stencil
!
 SUBROUTINE Trash_Stencil( myStencil )

   IMPLICIT NONE
   CLASS( Stencil ),INTENT(inout) :: myStencil

   DEALLOCATE( myStencil % relativeNeighbors )

 END SUBROUTINE Trash_Stencil

END MODULE POP_Stencil_Class
