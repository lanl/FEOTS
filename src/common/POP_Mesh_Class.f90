! POP_Mesh_Class.f90
! 
! Copyright 2016 Joseph Schoonover, Los Alamos National Laboratory (jschoonover@lanl.gov) 
! All rights reserved. 
! 
! POP_Mesh_Class.f90 is part of the Fast Equilibration of Ocean Tracers Software (FEOTS). 
! 
! //////////////////////////////////////////////////////////////////////////////////////////////// !
 
 

MODULE POP_Mesh_Class
!
!  This module was written as part of an offline diagnostic for the Fast Equilibrium of Ocean 
!  Tracers Software (FEOTS). This module provides wet-point to single dimension array mapping that is
!  used in the various solver modules.
! ================================================================================================ !

USE ModelPrecision
USE ConstantsDictionary
USE BinaryIO
! External Libraries
USE netcdf 


 IMPLICIT NONE


   TYPE POP_Mesh
      INTEGER                 :: nX, nY, nZ, nDOF
      INTEGER                 :: meshType
      ! Tracer mesh
      REAL(prec), ALLOCATABLE :: tLon(:,:), tLat(:,:)
      REAL(prec), ALLOCATABLE :: dXt(:,:), dYt(:,:)
      REAL(prec), ALLOCATABLE :: tArea(:,:)
      INTEGER, ALLOCATABLE    :: KmT(:,:)
      ! Vertical grid
      REAL(prec), ALLOCATABLE :: z(:), dz(:)
      ! Derived quantities
      REAL(prec), ALLOCATABLE :: tracermask(:,:,:)
      INTEGER, ALLOCATABLE    :: DOFtoIJK(:,:), IJKtoDOF(:,:,:)

       
       CONTAINS

          PROCEDURE :: Build => Build_POP_Mesh
          PROCEDURE :: Trash => Trash_POP_Mesh
          PROCEDURE :: Load  => Load_POP_Mesh
          
          PROCEDURE :: ConstructDummyTracerMesh => ConstructDummyTracerMesh_POP_Mesh
          PROCEDURE :: ConstructWetPointMap     => ConstructWetPointMap_POP_Mesh
          PROCEDURE :: MapFromDOFtoIJK          => MapFromDOFtoIJK_POP_Mesh
          PROCEDURE :: MapFromIJKtoDOF          => MapFromIJKtoDOF_POP_Mesh


    END TYPE POP_Mesh


 CONTAINS
!
!==================================================================================================!
!----------------------------- Manual Constructor/Destructor --------------------------------------!
!==================================================================================================!
!
 SUBROUTINE Build_POP_Mesh( theGrid, nX, nY, nZ )
 ! S/R Build
 !  
 !    
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(POP_Mesh), INTENT(out) :: theGrid
   INTEGER, INTENT(in)         :: nX, nY, nZ

      theGrid % nX = nX
      theGrid % nY = nY
      theGrid % nZ = nZ 
      ! Tracer mesh
      ALLOCATE( theGrid % tLon(1:nX,1:nY), theGrid % tLat(1:nX,1:nY) )
      ALLOCATE( theGrid % dXt(1:nX,1:nY), theGrid % dYt(1:nX,1:nY) )
      ALLOCATE( theGrid % tArea(1:nX,1:nY) )
      ALLOCATE( theGrid % KmT(1:nX,1:nY) )
      ! Vertical Grid
      ALLOCATE( theGrid % z(1:nZ), theGrid % dz(1:nZ) )
      ! Derived Quantities
      ALLOCATE( theGrid % tracermask(1:nX,1:nY,1:nZ), theGrid % IJKtoDOF(1:nX,1:nY,1:nZ) )

      theGrid % KmT = nZ
      
 END SUBROUTINE Build_POP_Mesh
!
 SUBROUTINE Trash_POP_Mesh( theGrid )
 ! S/R Trash
 !  
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(POP_Mesh), INTENT(inout) :: theGrid


      DEALLOCATE( theGrid % tLon, theGrid % tLat )
      DEALLOCATE( theGrid % dXt, theGrid % dYt )
      DEALLOCATE( theGrid % tArea )
      DEALLOCATE( theGrid % KmT )
      DEALLOCATE( theGrid % z, theGrid % dz )
      DEALLOCATE( theGrid % tracermask )
      DEALLOCATE( theGrid % DOFtoIJK, theGrid % IJKtoDOF )

 END SUBROUTINE Trash_POP_Mesh
!
 SUBROUTINE Load_POP_Mesh( theGrid, ncFile )
 ! S/R Load
 !  Desription:
 !  
 !    Loads in the MITgcm grid where all of the standard grid files are
 !    located in dataDir. This subroutine must be called after the 
 !    data-sturcture has been built with S/R BuildPOP_Mesh
 !    
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(POP_Mesh), INTENT(inout) :: theGrid
   CHARACTER(*), INTENT(in)       :: ncFile
   ! Local
   INTEGER       :: ncstatus, ncid, varid, dimid
   INTEGER       :: nX, nY, nZ
   CHARACTER(10) :: dimname

        ! ** Need to switch this to NetCDF ** 
        PRINT*, 'S/R Load_POP_Mesh : Reading in the grid information from '//TRIM(ncFile)
        
        ! Open the netcdf file - store the file handle in ncid        
        ncstatus = nf90_open( TRIM(ncFile), NF90_NOWRITE, ncid ) 
        CALL Check( ncstatus )
        ! Get the dimensions of the mesh !
        ncstatus = nf90_inq_dimid( ncid, "nlon", dimid )
        CALL Check( ncstatus )
        ncstatus = nf90_inquire_dimension( ncid, dimid, dimname, nX )
        CALL Check( ncstatus )
        
        ncstatus = nf90_inq_dimid( ncid, "nlat", dimid )
        CALL Check( ncstatus )
        ncstatus = nf90_inquire_dimension( ncid, dimid, dimname, nY )
        CALL Check( ncstatus )
        
        ncstatus = nf90_inq_dimid( ncid, "z_t", dimid )
        CALL Check( ncstatus )
        ncstatus = nf90_inquire_dimension( ncid, dimid, dimname, nZ )
        CALL Check( ncstatus )

        PRINT*, 'S/R Load_POP_Mesh : Grid Dimensions (nX,nY,nZ) : (',nX,',',nY,',',nZ,')'
        CALL theGrid % Build( nX, nY, nZ )
        
        ! Get the variable ID for the longitude
        ncstatus = nf90_inq_varid( ncid, "TLONG", varid )
        CALL Check( ncstatus )
        ! And read the longitude data in
        ncstatus = nf90_get_var( ncid, varid, theGrid % tLon )
        CALL Check( ncstatus )
        ! Get the variable ID for the latitude
        ncstatus = nf90_inq_varid( ncid, "TLAT", varid )
        CALL Check( ncstatus )
        ! And read the longitude data in
        ncstatus = nf90_get_var( ncid, varid, theGrid % tLat )
        CALL Check( ncstatus )

        ! Get the variable ID for the KMT
        ncstatus = nf90_inq_varid( ncid, "KMT", varid )
        CALL Check( ncstatus )
        ! And read the longitude data in
        ncstatus = nf90_get_var( ncid, varid, theGrid % kmt )
        CALL Check( ncstatus )

        ! Get the variable ID for the vertical grid
        ncstatus = nf90_inq_varid( ncid, "z_t", varid )
        CALL Check( ncstatus )
        ! And read the longitude data in
        ncstatus = nf90_get_var( ncid, varid, theGrid % z )
        CALL Check( ncstatus )

        ! Get the variable ID for the vertical grid spacing
        ncstatus = nf90_inq_varid( ncid, "dz", varid )
        CALL Check( ncstatus )
        ! And read the longitude data in
        ncstatus = nf90_get_var( ncid, varid, theGrid % dz )
        CALL Check( ncstatus )

        ! Get the variable ID for the lateral grid spacing in x
        ncstatus = nf90_inq_varid( ncid, "DXT", varid )
        CALL Check( ncstatus )
        ! And read the longitude data in
        ncstatus = nf90_get_var( ncid, varid, theGrid % dxt )
        CALL Check( ncstatus )

        ! Get the variable ID for the lateral grid spacing in y
        ncstatus = nf90_inq_varid( ncid, "DYT", varid )
        CALL Check( ncstatus )
        ! And read the longitude data in
        ncstatus = nf90_get_var( ncid, varid, theGrid % dyt )
        CALL Check( ncstatus )

        ! Get the variable ID for the lateral grid cell areas
        ncstatus = nf90_inq_varid( ncid, "TAREA", varid )
        CALL Check( ncstatus )
        ! And read the longitude data in
        ncstatus = nf90_get_var( ncid, varid, theGrid % tArea )
        CALL Check( ncstatus )

        ! Close the netcdf file
        ncstatus = nf90_close( ncid )
        CALL Check( ncstatus )


        CALL theGrid % ConstructWetPointMap( )
        
 END SUBROUTINE Load_POP_Mesh
!
!
!
 SUBROUTINE ConstructDummyTracerMesh_POP_Mesh( theGrid )
 ! S/R ConstructDummyTracerMesh
 !
 ! Builds a tracer mesh between (0,1)X(0,1)X(-1,0).
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(POP_Mesh), INTENT(inout) :: theGrid
   ! Local
   INTEGER    :: iX, iY, iZ
   REAL(prec) :: dx, dy, dz


      dx = ONE/REAL(theGrid % nX,prec)
      dy = ONE/REAL(theGrid % nY,prec)
      dz = ONE/REAL(theGrid % nZ,prec)

      theGrid % dxT = dx
      theGrid % dyT = dy
      theGrid % dz  = dz
      theGrid % tArea = dx*dy

      DO iX = 1, theGrid % nX
         theGrid % tLon(iX,:) = dx*( REAL(iX-1,prec) + HALF ) 
      ENDDO

      DO iY = 1, theGrid % nY
         theGrid % tLat(:,iY) = dy*( REAL(iY-1,prec) + HALF ) 
      ENDDO

      DO iZ = 1, theGrid % nZ
         theGrid % z(iZ) = dz*( REAL(iZ-1,prec) + HALF ) - ONE 
      ENDDO

 END SUBROUTINE ConstructDummyTracerMesh_POP_Mesh
!
!==================================================================================================!
!-------------------------------------- Type Specific ---------------------------------------------!
!==================================================================================================!
!
 SUBROUTINE ConstructWetPointMap_POP_Mesh( theGrid )
 ! S/R ConstructWetPointMap
 !
 !  This subroutine uses the KMT field to determine the wet-points. First, the "tracermask" is 
 !  set to one for all of the wet-points, and is left as zero for "dry" points in the structured
 !  "ijk" POP mesh. The sum of the tracermask gives the number of wet-points (aka "degrees of 
 !  freedom", or "DOF" ). Then, the tracermask is used to generate an array that is indexed over
 !  the DOF and returns the associated (i,j,k) in the original POP mesh.
 !
 !  ** This routine is used by "Build" and assumes that memory has already been allocated for 
 !     the POPmesh attributes. Since this routine is called by Build and should really not be 
 !     called again, this subroutine is PRIVATE ** 
 !  
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( POP_Mesh ), INTENT(inout) :: theGrid
   ! LOCAL
   INTEGER :: i, j, k, nDOF

      nDOF = 0
      theGrid % tracerMask = ZERO

      ! Count the number of wet-points and set the IJK to DOF mapping
      DO j = 1, theGrid % nY
         DO i = 1, theGrid % nX

            DO k = 1, theGrid % KMT(i,j)
               nDOF = nDOF + 1
               theGrid % tracerMask(i,j,k) = ONE
               theGrid % IJKtoDOF(i,j,k)   = nDOF
            ENDDO

         ENDDO
      ENDDO

      PRINT*, 'S/R ConstructWetPointMap :'
      PRINT*, 'Found ', nDOF, 'degrees of freedom from ', theGrid % nX*theGrid % nY*theGrid % nZ, 'mesh points.'

      ALLOCATE( theGrid % DOFtoIJK(1:3,1:nDOF) ) 
      theGrid % nDOF = nDOF

      ! Now we can set the DOF to IJK mapping
      nDOF = 0
      DO j = 1, theGrid % nY
         DO i = 1, theGrid % nX

            DO k = 1, theGrid % KMT(i,j)
               nDOF = nDOF + 1
               theGrid % DOFtoIJK(1,nDOF) = i
               theGrid % DOFtoIJK(2,nDOF) = j
               theGrid % DOFtoIJK(3,nDOF) = k
            ENDDO

         ENDDO
      ENDDO

 END SUBROUTINE ConstructWetPointMap_POP_Mesh
!
!
! 
 FUNCTION MapFromDOFtoIJK_POP_Mesh( theGrid, dofArray ) RESULT( ijkArray )
 ! Function MapFromDOFtoIJK
 !
 !   This function takes the single-dimension array, indexed from 1 to nDOF, and maps it back to
 !   the structured POP storage, a 3-D array. The POP_Mesh attribute "DOFtoIJK" facilitates this
 !   mapping and is assumed to already be constructed before calling this function. If a call
 !   to "Build" has been issued for "theGrid", then the mapping has already been created.
 !
 !   Dry points are filled in with the specified "fvalue"
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( POP_Mesh ) :: theGrid
   REAL(prec)       :: dofArray(1:theGrid % nDOF)
   REAL(prec)       :: ijkArray(1:theGrid % nX, 1:theGrid % nY, 1:theGrid % nZ)
   ! Local
   INTEGER :: i, j, k, l


      ijkArray = fillValue
      DO l = 1, theGrid % nDof

         i = theGrid % DOFtoIJK(l,1)
         j = theGrid % DOFtoIJK(l,2)
         k = theGrid % DOFtoIJK(l,3)

         ijkArray(i,j,k) = dofArray(l)

      ENDDO

 END FUNCTION MapFromDOFtoIJK_POP_Mesh
!
!
!
 FUNCTION MapFromIJKtoDOF_POP_Mesh( theGrid, ijkArray ) RESULT( dofArray )
 ! Function MapFromIJKtoDOF
 !
 !   This function takes the 3-D "POP-storage" array and maps it to a single dimensioned "DOF" array.
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( POP_Mesh ) :: theGrid
   REAL(prec)       :: ijkArray(1:theGrid % nX, 1:theGrid % nY, 1:theGrid % nZ)
   REAL(prec)       :: dofArray(1:theGrid % nDOF)
   ! Local
   INTEGER :: i, j, k, l

      ijkArray = fillValue
      DO l = 1, theGrid % nDof

         i = theGrid % DOFtoIJK(l,1)
         j = theGrid % DOFtoIJK(l,2)
         k = theGrid % DOFtoIJK(l,3)

         dofArray(l) = ijkArray(i,j,k)

      ENDDO

 END FUNCTION MapFromIJKtoDOF_POP_Mesh
!
!
!
! SUBROUTINE GenerateRegionalMasks_POP_Mesh( theGrid )
   ! Generates the (i,j) and associated DOF indices that fall within a latitude and longitude
   ! bounded box. This provides necessary information for performing regional simulations.
   !
! END SUBROUTINE GenerateRegionalMasks_POP_Mesh
!
!==================================================================================================!
!---------------------------------------- File I/O ------------------------------------------------!
!==================================================================================================!
!
  SUBROUTINE Check(status)
    INTEGER, INTENT (in) :: status
    
    IF(status /= nf90_noerr) THEN 
      PRINT *, trim(nf90_strerror(status))
      STOP "NetCDF Error, Stopped"
    ENDIF
  END SUBROUTINE Check 
!
!
!
END MODULE POP_Mesh_Class
