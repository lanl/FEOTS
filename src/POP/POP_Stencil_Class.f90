! POP_Stencil_Class.f90
! 
! Copyright 2016 Joseph Schoonover, Los Alamos National Laboratory (jschoonover@lanl.gov) 
! All rights reserved. 
! 
! POP_Stencil_Class.f90 is part of the Fast Equilibration of Ocean Tracers Software (FEOTS). 
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !
 
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

      IF( stencilFlag == LaxWendroff .AND. flavor == Overlap )THEN
         
         PRINT*, 'S/R : Build_Stencil : Constructing Lax-Wendroff with overlap flavor.'
         myStencil % nPoints = 24
         ALLOCATE( myStencil % relativeNeighbors(1:3,1:24) )

         ! 
         myStencil % relativeNeighbors(1,1) = 0 ! x
         myStencil % relativeNeighbors(2,1) = -2  ! y
         myStencil % relativeNeighbors(3,1) = 0  ! z
         !
         myStencil % relativeNeighbors(1,2) = -1  ! x
         myStencil % relativeNeighbors(2,2) = -1  ! y
         myStencil % relativeNeighbors(3,2) = -0  ! z
         ! 
         myStencil % relativeNeighbors(1,3) = 0  ! x
         myStencil % relativeNeighbors(2,3) = -1  ! y
         myStencil % relativeNeighbors(3,3) = 0  ! z
         ! 
         myStencil % relativeNeighbors(1,4) = 1  ! x
         myStencil % relativeNeighbors(2,4) = -1 ! y
         myStencil % relativeNeighbors(3,4) = 0  ! z
         ! 
         myStencil % relativeNeighbors(1,5) = -2  ! x
         myStencil % relativeNeighbors(2,5) = 0  ! y
         myStencil % relativeNeighbors(3,5) = 0  ! z
         ! 
         myStencil % relativeNeighbors(1,6) = -1  ! x
         myStencil % relativeNeighbors(2,6) = 0  ! y
         myStencil % relativeNeighbors(3,6) = 0 ! z
         ! 
         myStencil % relativeNeighbors(1,7) = 1  ! x
         myStencil % relativeNeighbors(2,7) = 0  ! y
         myStencil % relativeNeighbors(3,7) = 0  ! z

         myStencil % relativeNeighbors(1,8) = 2 ! x
         myStencil % relativeNeighbors(2,8) = 0  ! y
         myStencil % relativeNeighbors(3,8) = 0  ! z
         
         myStencil % relativeNeighbors(1,9) = -1 ! x
         myStencil % relativeNeighbors(2,9) = 1  ! y
         myStencil % relativeNeighbors(3,9) = 0  ! z
         
         myStencil % relativeNeighbors(1,10) = 0 ! x
         myStencil % relativeNeighbors(2,10) = 1  ! y
         myStencil % relativeNeighbors(3,10) = 0  ! z
         
         myStencil % relativeNeighbors(1,11) = 1 ! x
         myStencil % relativeNeighbors(2,11) = 1  ! y
         myStencil % relativeNeighbors(3,11) = 0  ! z
         
         myStencil % relativeNeighbors(1,12) = 0 ! x
         myStencil % relativeNeighbors(2,12) = 2  ! y
         myStencil % relativeNeighbors(3,12) = 0  ! z
         
         myStencil % relativeNeighbors(1,13) = 0 ! x
         myStencil % relativeNeighbors(2,13) = -1  ! y
         myStencil % relativeNeighbors(3,13) = -1  ! z
         
         myStencil % relativeNeighbors(1,14) = -1 ! x
         myStencil % relativeNeighbors(2,14) = 0  ! y
         myStencil % relativeNeighbors(3,14) = -1  ! z

         myStencil % relativeNeighbors(1,15) = 1 ! x
         myStencil % relativeNeighbors(2,15) = 0  ! y
         myStencil % relativeNeighbors(3,15) = -1  ! z
         
         myStencil % relativeNeighbors(1,16) = 0 ! x
         myStencil % relativeNeighbors(2,16) = 1  ! y
         myStencil % relativeNeighbors(3,16) = -1  ! z
         
         myStencil % relativeNeighbors(1,17) = 0 ! x
         myStencil % relativeNeighbors(2,17) = 0  ! y
         myStencil % relativeNeighbors(3,17) = -1  ! z

         myStencil % relativeNeighbors(1,18) = 0  ! x
         myStencil % relativeNeighbors(2,18) = 0  ! y
         myStencil % relativeNeighbors(3,18) = -2  ! z
         ! Below -- 1 level
         myStencil % relativeNeighbors(1,19) = 0  ! x
         myStencil % relativeNeighbors(2,19) = -1  ! y
         myStencil % relativeNeighbors(3,19) = 1  ! z
         ! 
         myStencil % relativeNeighbors(1,20) = -1  ! x
         myStencil % relativeNeighbors(2,20) = 0 ! y
         myStencil % relativeNeighbors(3,20) = 1  ! z
         ! 
         myStencil % relativeNeighbors(1,21) = 1  ! x
         myStencil % relativeNeighbors(2,21) = 0 ! y
         myStencil % relativeNeighbors(3,21) = 1  ! z
         ! 
         myStencil % relativeNeighbors(1,22) = 0  ! x
         myStencil % relativeNeighbors(2,22) = 1  ! y
         myStencil % relativeNeighbors(3,22) = 1  ! z
         ! 
         myStencil % relativeNeighbors(1,23) = 0  ! x
         myStencil % relativeNeighbors(2,23) = 0  ! y
         myStencil % relativeNeighbors(3,23) = 1 ! z
         ! Below -- 2 levels
         myStencil % relativeNeighbors(1,24) = 0 ! x
         myStencil % relativeNeighbors(2,24) = 0  ! y
         myStencil % relativeNeighbors(3,24) = 2  ! z

!         ! Above -- 2 levels
!         myStencil % relativeNeighbors(1,1) = 0 ! x
!         myStencil % relativeNeighbors(2,1) = 0  ! y
!         myStencil % relativeNeighbors(3,1) = -2  ! z
!         ! -----------------
!         ! Above -- 1 level
!         myStencil % relativeNeighbors(1,2) = 0  ! x
!         myStencil % relativeNeighbors(2,2) = -1  ! y
!         myStencil % relativeNeighbors(3,2) = -1  ! z
!         ! 
!         myStencil % relativeNeighbors(1,3) = -1  ! x
!         myStencil % relativeNeighbors(2,3) = 0  ! y
!         myStencil % relativeNeighbors(3,3) = -1  ! z
!         ! 
!         myStencil % relativeNeighbors(1,4) = 0  ! x
!         myStencil % relativeNeighbors(2,4) = 0 ! y
!         myStencil % relativeNeighbors(3,4) = -1  ! z
!         ! 
!         myStencil % relativeNeighbors(1,5) = 1  ! x
!         myStencil % relativeNeighbors(2,5) = 0  ! y
!         myStencil % relativeNeighbors(3,5) = -1  ! z
!         ! 
!         myStencil % relativeNeighbors(1,6) = 0  ! x
!         myStencil % relativeNeighbors(2,6) = 1  ! y
!         myStencil % relativeNeighbors(3,6) = -1 ! z
!         ! -----------------
!         ! At level
!         myStencil % relativeNeighbors(1,7) = 0  ! x
!         myStencil % relativeNeighbors(2,7) = -2  ! y
!         myStencil % relativeNeighbors(3,7) = 0  ! z
!
!         myStencil % relativeNeighbors(1,8) = -1 ! x
!         myStencil % relativeNeighbors(2,8) = -1  ! y
!         myStencil % relativeNeighbors(3,8) = 0  ! z
!         
!         myStencil % relativeNeighbors(1,9) = 0 ! x
!         myStencil % relativeNeighbors(2,9) = -1  ! y
!         myStencil % relativeNeighbors(3,9) = 0  ! z
!         
!         myStencil % relativeNeighbors(1,10) = 1 ! x
!         myStencil % relativeNeighbors(2,10) = -1  ! y
!         myStencil % relativeNeighbors(3,10) = 0  ! z
!         
!         myStencil % relativeNeighbors(1,11) = -2 ! x
!         myStencil % relativeNeighbors(2,11) = 0  ! y
!         myStencil % relativeNeighbors(3,11) = 0  ! z
!         
!         myStencil % relativeNeighbors(1,12) = -1 ! x
!         myStencil % relativeNeighbors(2,12) = 0  ! y
!         myStencil % relativeNeighbors(3,12) = 0  ! z
!         
!         myStencil % relativeNeighbors(1,13) = 1 ! x
!         myStencil % relativeNeighbors(2,13) = 0  ! y
!         myStencil % relativeNeighbors(3,13) = 0  ! z
!         
!         myStencil % relativeNeighbors(1,14) = 2 ! x
!         myStencil % relativeNeighbors(2,14) = 0  ! y
!         myStencil % relativeNeighbors(3,14) = 0  ! z
!
!         myStencil % relativeNeighbors(1,15) = -1 ! x
!         myStencil % relativeNeighbors(2,15) = 1  ! y
!         myStencil % relativeNeighbors(3,15) = 0  ! z
!         
!         myStencil % relativeNeighbors(1,16) = 0 ! x
!         myStencil % relativeNeighbors(2,16) = 1  ! y
!         myStencil % relativeNeighbors(3,16) = 0  ! z
!         
!         myStencil % relativeNeighbors(1,17) = 1 ! x
!         myStencil % relativeNeighbors(2,17) = 1  ! y
!         myStencil % relativeNeighbors(3,17) = 0  ! z
!
!         myStencil % relativeNeighbors(1,18) = 0  ! x
!         myStencil % relativeNeighbors(2,18) = 2  ! y
!         myStencil % relativeNeighbors(3,18) = 0  ! z
!         ! Below -- 1 level
!         myStencil % relativeNeighbors(1,19) = 0  ! x
!         myStencil % relativeNeighbors(2,19) = -1  ! y
!         myStencil % relativeNeighbors(3,19) = 1  ! z
!         ! 
!         myStencil % relativeNeighbors(1,20) = -1  ! x
!         myStencil % relativeNeighbors(2,20) = 0  ! y
!         myStencil % relativeNeighbors(3,20) = 1  ! z
!         ! 
!         myStencil % relativeNeighbors(1,21) = 0  ! x
!         myStencil % relativeNeighbors(2,21) = 0 ! y
!         myStencil % relativeNeighbors(3,21) = 1  ! z
!         ! 
!         myStencil % relativeNeighbors(1,22) = 1  ! x
!         myStencil % relativeNeighbors(2,22) = 0  ! y
!         myStencil % relativeNeighbors(3,22) = 1  ! z
!         ! 
!         myStencil % relativeNeighbors(1,23) = 0  ! x
!         myStencil % relativeNeighbors(2,23) = 1  ! y
!         myStencil % relativeNeighbors(3,23) = 1 ! z
!         ! Below -- 2 levels
!         myStencil % relativeNeighbors(1,24) = 0 ! x
!         myStencil % relativeNeighbors(2,24) = 0  ! y
!         myStencil % relativeNeighbors(3,24) = 2  ! z
          
      ELSEIF( stencilFlag == LaxWendroff .AND. flavor == Normal )THEN
         
         PRINT*, 'S/R : Build_Stencil : Constructing Lax-Wendroff stencil with Normal flavor.'
         myStencil % nPoints = 7
         ALLOCATE( myStencil % relativeNeighbors(1:3,1:7) )
         !
         myStencil % relativeNeighbors(1,1) = 0 ! x
         myStencil % relativeNeighbors(2,1) = 0  ! y
         myStencil % relativeNeighbors(3,1) = -1  ! z
         !
         myStencil % relativeNeighbors(1,2) = 0  ! x
         myStencil % relativeNeighbors(2,2) = 0  ! y
         myStencil % relativeNeighbors(3,2) = 1  ! z
         !
         myStencil % relativeNeighbors(1,3) = -1  ! x
         myStencil % relativeNeighbors(2,3) = 0  ! y
         myStencil % relativeNeighbors(3,3) = 0  ! z
         !
         myStencil % relativeNeighbors(1,4) = 1  ! x
         myStencil % relativeNeighbors(2,4) = 0 ! y
         myStencil % relativeNeighbors(3,4) = 0  ! z
         !
         myStencil % relativeNeighbors(1,5) = 0  ! x
         myStencil % relativeNeighbors(2,5) = -1  ! y
         myStencil % relativeNeighbors(3,5) = 0  ! z
         !
         myStencil % relativeNeighbors(1,6) = 0  ! x
         myStencil % relativeNeighbors(2,6) = 1  ! y
         myStencil % relativeNeighbors(3,6) = 0 ! z
         !
         myStencil % relativeNeighbors(1,7) = 0  ! x
         myStencil % relativeNeighbors(2,7) = 0  ! y
         myStencil % relativeNeighbors(3,7) = 0  ! z
!         ! west
!         myStencil % relativeNeighbors(1,1) = -1 ! x
!         myStencil % relativeNeighbors(2,1) = 0  ! y
!         myStencil % relativeNeighbors(3,1) = 0  ! z
!         ! central
!         myStencil % relativeNeighbors(1,2) = 0  ! x
!         myStencil % relativeNeighbors(2,2) = 0  ! y
!         myStencil % relativeNeighbors(3,2) = 0  ! z
!         ! east
!         myStencil % relativeNeighbors(1,3) = 1  ! x
!         myStencil % relativeNeighbors(2,3) = 0  ! y
!         myStencil % relativeNeighbors(3,3) = 0  ! z
!         ! south
!         myStencil % relativeNeighbors(1,4) = 0  ! x
!         myStencil % relativeNeighbors(2,4) = -1 ! y
!         myStencil % relativeNeighbors(3,4) = 0  ! z
!         ! north
!         myStencil % relativeNeighbors(1,5) = 0  ! x
!         myStencil % relativeNeighbors(2,5) = 1  ! y
!         myStencil % relativeNeighbors(3,5) = 0  ! z
!         ! above
!         myStencil % relativeNeighbors(1,6) = 0  ! x
!         myStencil % relativeNeighbors(2,6) = 0  ! y
!         myStencil % relativeNeighbors(3,6) = -1 ! z
!         ! below
!         myStencil % relativeNeighbors(1,7) = 0  ! x
!         myStencil % relativeNeighbors(2,7) = 0  ! y
!         myStencil % relativeNeighbors(3,7) = 1  ! z

      ELSEIF( stencilFlag == LaxWendroff .AND. flavor == Lateral )THEN
         
         PRINT*, 'S/R : Build_Stencil : Constructing Lax-Wendroff stencil with Lateral flavor.'
         myStencil % nPoints = 5
         ALLOCATE( myStencil % relativeNeighbors(1:3,1:5) )
         ! west
         myStencil % relativeNeighbors(1,1) = -1 ! x
         myStencil % relativeNeighbors(2,1) = 0  ! y
         myStencil % relativeNeighbors(3,1) = 0  ! z
         ! central
         myStencil % relativeNeighbors(1,2) = 0  ! x
         myStencil % relativeNeighbors(2,2) = 0  ! y
         myStencil % relativeNeighbors(3,2) = 0  ! z
         ! east
         myStencil % relativeNeighbors(1,3) = 1  ! x
         myStencil % relativeNeighbors(2,3) = 0  ! y
         myStencil % relativeNeighbors(3,3) = 0  ! z
         ! south
         myStencil % relativeNeighbors(1,4) = 0  ! x
         myStencil % relativeNeighbors(2,4) = -1 ! y
         myStencil % relativeNeighbors(3,4) = 0  ! z
         ! north
         myStencil % relativeNeighbors(1,5) = 0  ! x
         myStencil % relativeNeighbors(2,5) = 1  ! y
         myStencil % relativeNeighbors(3,5) = 0  ! z

      ELSEIF( stencilFlag == LaxWendroff .AND. flavor == LateralPlusCorners )THEN
         
         PRINT*, 'S/R : Build_Stencil : Constructing Lax-Wendroff stencil with LateralPlusCorners flavor.'
         myStencil % nPoints = 9
         ALLOCATE( myStencil % relativeNeighbors(1:3,1:9) )
         ! south-west
         myStencil % relativeNeighbors(1,1) = -1 ! x
         myStencil % relativeNeighbors(2,1) = -1  ! y
         myStencil % relativeNeighbors(3,1) = 0  ! z
         ! south
         myStencil % relativeNeighbors(1,2) = 0  ! x
         myStencil % relativeNeighbors(2,2) = -1 ! y
         myStencil % relativeNeighbors(3,2) = 0  ! z
         ! south-east
         myStencil % relativeNeighbors(1,3) = 1  ! x
         myStencil % relativeNeighbors(2,3) = -1 ! y
         myStencil % relativeNeighbors(3,3) = 0  ! z

         ! west
         myStencil % relativeNeighbors(1,4) = -1 ! x
         myStencil % relativeNeighbors(2,4) = 0  ! y
         myStencil % relativeNeighbors(3,4) = 0  ! z
         ! central
         myStencil % relativeNeighbors(1,5) = 0  ! x
         myStencil % relativeNeighbors(2,5) = 0 ! y
         myStencil % relativeNeighbors(3,5) = 0  ! z
         ! east
         myStencil % relativeNeighbors(1,6) = 1  ! x
         myStencil % relativeNeighbors(2,6) = 0 ! y
         myStencil % relativeNeighbors(3,6) = 0  ! z

         ! north-west
         myStencil % relativeNeighbors(1,7) = -1 ! x
         myStencil % relativeNeighbors(2,7) = 1  ! y
         myStencil % relativeNeighbors(3,7) = 0  ! z
         ! north
         myStencil % relativeNeighbors(1,8) = 0  ! x
         myStencil % relativeNeighbors(2,8) = 1 ! y
         myStencil % relativeNeighbors(3,8) = 0  ! z
         ! north-east
         myStencil % relativeNeighbors(1,9) = 1  ! x
         myStencil % relativeNeighbors(2,9) = 1 ! y
         myStencil % relativeNeighbors(3,9) = 0  ! z



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
