! POP_Mesh_Class.f90
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
 
 

MODULE POP_Mesh_Class
!
!  This module was written as part of an offline diagnostic for the Fast Equilibrium of Ocean 
!  Tracers Software (FEOTS). This module provides wet-point to single dimension array mapping that is
!  used in the various solver modules.
! ================================================================================================ !

USE ModelPrecision
USE ConstantsDictionary
USE BinaryIO
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
      REAL(prec), ALLOCATABLE :: z(:), dz(:), dzw(:)
      ! Derived quantities
      REAL(prec), ALLOCATABLE :: tracermask(:,:,:)
      INTEGER, ALLOCATABLE    :: DOFtoIJK(:,:), IJKtoDOF(:,:,:)

       
       CONTAINS

          PROCEDURE :: Build => Build_POP_Mesh
          PROCEDURE :: Trash => Trash_POP_Mesh
          PROCEDURE :: Load  => Load_POP_Mesh
          PROCEDURE :: LoadWithMask  => LoadWithMask_POP_Mesh
    
          PROCEDURE :: WriteNetCDF => WriteNetCDF_POP_Mesh          

          PROCEDURE :: ConstructDummyTracerMesh => ConstructDummyTracerMesh_POP_Mesh
          PROCEDURE :: ConstructWetPointMap     => ConstructWetPointMap_POP_Mesh
          PROCEDURE :: ConstructWetPointMapWithMask     => ConstructWetPointMapWithMask_POP_Mesh
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
      ALLOCATE( theGrid % z(1:nZ), theGrid % dz(1:nZ), theGrid % dzw(1:nZ) )
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
      DEALLOCATE( theGrid % z, theGrid % dz, theGrid % dzw )
      DEALLOCATE( theGrid % tracermask )
      IF( ALLOCATED( theGrid % DOFtoIJK ) )THEN
         DEALLOCATE( theGrid % DOFtoIJK, theGrid % IJKtoDOF )
      ENDIF
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
   INTEGER       :: ncid, varid, dimid, nstat
   INTEGER       :: nX, nY, nZ
   CHARACTER(10) :: dimname

        ! ** Need to switch this to NetCDF ** 
        PRINT*, 'S/R Load_POP_Mesh : Reading in the grid information from '//TRIM(ncFile)
        
        ! Open the netcdf file - store the file handle in ncid        
        CALL Check( nf90_open( TRIM(ncFile), NF90_NOWRITE, ncid ) ) 
        ! Get the dimensions of the mesh !
        CALL Check( nf90_inq_dimid( ncid, "nlon", dimid ) ) 
        CALL Check( nf90_inquire_dimension( ncid, dimid, dimname, nX ) )
        
        CALL Check( nf90_inq_dimid( ncid, "nlat", dimid ) )
        CALL Check( nf90_inquire_dimension( ncid, dimid, dimname, nY ) )
        
        CALL Check( nf90_inq_dimid( ncid, "z_t", dimid ) ) 
        CALL Check( nf90_inquire_dimension( ncid, dimid, dimname, nZ ) )

        PRINT*, 'S/R Load_POP_Mesh : Grid Dimensions (nX,nY,nZ) : (',nX,',',nY,',',nZ,')'
        CALL theGrid % Build( nX, nY, nZ )
        
        ! Get the variable ID for the longitude
        CALL Check( nf90_inq_varid( ncid, "TLONG", varid ) )
        ! And read the longitude data in
        CALL Check( nf90_get_var( ncid, varid, theGrid % tLon ) )
        ! Get the variable ID for the latitude
        CALL Check( nf90_inq_varid( ncid, "TLAT", varid ) )
        ! And read the longitude data in
        CALL Check( nf90_get_var( ncid, varid, theGrid % tLat ) )

        ! Get the variable ID for the KMT
        CALL Check( nf90_inq_varid( ncid, "KMT", varid ) )
        ! And read the longitude data in
        CALL Check( nf90_get_var( ncid, varid, theGrid % kmt ) )

        ! Get the variable ID for the vertical grid
        CALL Check( nf90_inq_varid( ncid, "z_t", varid ) )
        ! And read the longitude data in
        CALL Check( nf90_get_var( ncid, varid, theGrid % z ) )

        ! Get the variable ID for the vertical grid spacing
        CALL Check( nf90_inq_varid( ncid, "dz", varid ) )
        ! And read the longitude data in
        CALL Check( nf90_get_var( ncid, varid, theGrid % dz ) )

        ! Get the variable ID for the vertical grid spacing
        CALL Check( nf90_inq_varid( ncid, "dzw", varid ) )
        ! And read the longitude data in
        CALL Check( nf90_get_var( ncid, varid, theGrid % dzw ) )


        ! Get the variable ID for the lateral grid spacing in x
        CALL Check( nf90_inq_varid( ncid, "DXT", varid ) )
        ! And read the longitude data in
        CALL Check( nf90_get_var( ncid, varid, theGrid % dxt ) )

        ! Get the variable ID for the lateral grid spacing in y
        CALL Check( nf90_inq_varid( ncid, "DYT", varid ) )
        ! And read the longitude data in
        CALL Check( nf90_get_var( ncid, varid, theGrid % dyt ) )

        ! Get the variable ID for the lateral grid cell areas
        CALL Check( nf90_inq_varid( ncid, "TAREA", varid ) )
        ! And read the longitude data in
        CALL Check( nf90_get_var( ncid, varid, theGrid % tArea ) )

        ! Get the variable ID for the lateral grid cell areas
        nstat = nf90_inq_varid( ncid, "TracerMask", varid )
        IF( nstat == nf90_noerr )THEN
           CALL Check( nf90_get_var( ncid, varid, theGrid % tracermask ) )
        ENDIF
        ! Close the netcdf file
        CALL Check( nf90_close( ncid ) )

        CALL theGrid % ConstructWetPointMap( )
        
 END SUBROUTINE Load_POP_Mesh
!
 SUBROUTINE LoadWithMask_POP_Mesh( theGrid, ncFile )
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
   INTEGER       :: ncid, varid, dimid, nstat
   INTEGER       :: nX, nY, nZ
   CHARACTER(10) :: dimname

        ! ** Need to switch this to NetCDF ** 
        PRINT*, 'S/R Load_POP_Mesh : Reading in the grid information from '//TRIM(ncFile)
        
        ! Open the netcdf file - store the file handle in ncid        
        CALL Check( nf90_open( TRIM(ncFile), NF90_NOWRITE, ncid ) ) 
        ! Get the dimensions of the mesh !
        CALL Check( nf90_inq_dimid( ncid, "nlon", dimid ) ) 
        CALL Check( nf90_inquire_dimension( ncid, dimid, dimname, nX ) )
        
        CALL Check( nf90_inq_dimid( ncid, "nlat", dimid ) )
        CALL Check( nf90_inquire_dimension( ncid, dimid, dimname, nY ) )
        
        CALL Check( nf90_inq_dimid( ncid, "z_t", dimid ) ) 
        CALL Check( nf90_inquire_dimension( ncid, dimid, dimname, nZ ) )

        PRINT*, 'S/R Load_POP_Mesh : Grid Dimensions (nX,nY,nZ) : (',nX,',',nY,',',nZ,')'
        CALL theGrid % Build( nX, nY, nZ )
        
        ! Get the variable ID for the longitude
        CALL Check( nf90_inq_varid( ncid, "TLONG", varid ) )
        ! And read the longitude data in
        CALL Check( nf90_get_var( ncid, varid, theGrid % tLon ) )
        ! Get the variable ID for the latitude
        CALL Check( nf90_inq_varid( ncid, "TLAT", varid ) )
        ! And read the longitude data in
        CALL Check( nf90_get_var( ncid, varid, theGrid % tLat ) )

        ! Get the variable ID for the KMT
        CALL Check( nf90_inq_varid( ncid, "KMT", varid ) )
        ! And read the longitude data in
        CALL Check( nf90_get_var( ncid, varid, theGrid % kmt ) )

        ! Get the variable ID for the vertical grid
        CALL Check( nf90_inq_varid( ncid, "z_t", varid ) )
        ! And read the longitude data in
        CALL Check( nf90_get_var( ncid, varid, theGrid % z ) )

        ! Get the variable ID for the vertical grid spacing
        CALL Check( nf90_inq_varid( ncid, "dz", varid ) )
        ! And read the longitude data in
        CALL Check( nf90_get_var( ncid, varid, theGrid % dz ) )

        ! Get the variable ID for the vertical grid spacing
        CALL Check( nf90_inq_varid( ncid, "dzw", varid ) )
        ! And read the longitude data in
        CALL Check( nf90_get_var( ncid, varid, theGrid % dzw ) )


        ! Get the variable ID for the lateral grid spacing in x
        CALL Check( nf90_inq_varid( ncid, "DXT", varid ) )
        ! And read the longitude data in
        CALL Check( nf90_get_var( ncid, varid, theGrid % dxt ) )

        ! Get the variable ID for the lateral grid spacing in y
        CALL Check( nf90_inq_varid( ncid, "DYT", varid ) )
        ! And read the longitude data in
        CALL Check( nf90_get_var( ncid, varid, theGrid % dyt ) )

        ! Get the variable ID for the lateral grid cell areas
        CALL Check( nf90_inq_varid( ncid, "TAREA", varid ) )
        ! And read the longitude data in
        CALL Check( nf90_get_var( ncid, varid, theGrid % tArea ) )

        ! Get the variable ID for the lateral grid cell areas
        nstat = nf90_inq_varid( ncid, "TracerMask", varid )
        IF( nstat == nf90_noerr )THEN
           CALL Check( nf90_get_var( ncid, varid, theGrid % tracermask ) )
        ENDIF
        ! Close the netcdf file
        CALL Check( nf90_close( ncid ) )

        CALL theGrid % ConstructWetPointMapWithMask( )
        
 END SUBROUTINE LoadWithMask_POP_Mesh
!
 SUBROUTINE WriteNetCDF_POP_Mesh( mesh, ncFile )
 ! S/R WriteNetCDF
   IMPLICIT NONE
   CLASS(POP_Mesh), INTENT(inout) :: mesh
   CHARACTER(*), INTENT(in)       :: ncFile
   ! Local
   INTEGER :: ncid, z_dimid, x_dimid, y_dimid
   INTEGER :: z_varid, x_varid, y_varid, kmt_varid
   INTEGER :: dz_varid, dzw_varid, dxt_varid, dyt_varid, tarea_varid, tmask_varid

        ! ** Need to switch this to NetCDF ** 
        PRINT*, 'S/R WriteNetCDF_POP_Mesh : Writing the grid information to '//TRIM(ncFile)
        
        ! Open the netcdf file - store the file handle in ncid        
        CALL Check( nf90_create( TRIM(ncFile), NF90_CLOBBER, ncid ) ) 
        PRINT*, '                           Defining dimensions of the mesh'
        ! Get the dimensions of the mesh !
        CALL Check( nf90_def_dim( ncid, "z_t", mesh % nZ, z_dimid ) ) 
        CALL Check( nf90_def_dim( ncid, "nlon", mesh % nX, x_dimid ) ) 
        CALL Check( nf90_def_dim( ncid, "nlat", mesh % nY, y_dimid ) ) 
       
        PRINT*,'                           Defining mesh variables' 
! Create variables -- here we need to create arrays for the dimensions
        CALL Check( nf90_def_var( ncid, "z_t", &
                                  NF90_DOUBLE, &
                                  z_dimid, z_varid ) )

        CALL Check( nf90_def_var( ncid, "TLAT", &
                                  NF90_DOUBLE, (/ x_dimid, y_dimid /),&
                                  y_varid ) )

        CALL Check( nf90_def_var( ncid, "TLONG", &
                                  NF90_DOUBLE, (/ x_dimid, y_dimid /),&
                                  x_varid ) )        

        CALL Check( nf90_def_var( ncid, "KMT", &
                                  NF90_INT, (/ x_dimid, y_dimid /),&
                                  kmt_varid ) )

        CALL Check( nf90_def_var( ncid, "dz", &
                                  NF90_DOUBLE, &
                                  z_dimid, dz_varid ) )

        CALL Check( nf90_def_var( ncid, "dzw", &
                                  NF90_DOUBLE, &
                                  z_dimid, dzw_varid ) )

        CALL Check( nf90_def_var( ncid, "DXT", &
                                  NF90_DOUBLE, (/ x_dimid, y_dimid /),&
                                  dxt_varid ) )        

        CALL Check( nf90_def_var( ncid, "DYT", &
                                  NF90_DOUBLE, (/ x_dimid, y_dimid /),&
                                  dyt_varid ) )        

        CALL Check( nf90_def_var( ncid, "TAREA", &
                                  NF90_DOUBLE, (/ x_dimid, y_dimid /),&
                                  tarea_varid ) )        
        
        CALL Check( nf90_def_var( ncid, "TracerMask", &
                                  NF90_DOUBLE, (/ x_dimid, y_dimid, z_dimid /),&
                                  tmask_varid ) )        
        
      CALL Check( nf90_put_att( ncid, tmask_varid, "long_name", "Domain Mask" ) )
      CALL Check( nf90_put_att( ncid, tmask_varid, "units", "" ) )
      CALL Check( nf90_put_att( ncid, tmask_varid, "_FillValue", fillValue) )
      CALL Check( nf90_put_att( ncid, tmask_varid, "missing_value", fillValue) )

        PRINT*,'                           Defining units.' 
        CALL Check( nf90_put_att( ncid, z_varid, 'units', 'cm' ) )
        CALL Check( nf90_put_att( ncid, y_varid, 'units', 'degrees North' ) )
        CALL Check( nf90_put_att( ncid, x_varid, 'units', 'degrees East' ) )
        CALL Check( nf90_put_att( ncid, kmt_varid, 'units', 'index' ) )
        CALL Check( nf90_put_att( ncid, dz_varid, 'units', 'cm' ) )
        CALL Check( nf90_put_att( ncid, dzw_varid, 'units', 'cm' ) )
        CALL Check( nf90_put_att( ncid, dxt_varid, 'units', 'km' ) )
        CALL Check( nf90_put_att( ncid, dyt_varid, 'units', 'km' ) )
        CALL Check( nf90_put_att( ncid, tarea_varid, 'units', 'L^2' ) )
        CALL Check( nf90_put_att( ncid, tmask_varid, 'units', 'unitless' ) )


        CALL Check( nf90_enddef(ncid) )
        PRINT*,'                           Writing variables to file.' 
        ! Put variables into the netcdf file
        CALL Check( nf90_put_var( ncid, z_varid, mesh % z ) )
        CALL Check( nf90_put_var( ncid, y_varid, mesh % tLat ) )
        CALL Check( nf90_put_var( ncid, x_varid, mesh % tLon ) )
        CALL Check( nf90_put_var( ncid, kmt_varid, mesh % kmt ) )
        CALL Check( nf90_put_var( ncid, dz_varid, mesh % dz ) )
        CALL Check( nf90_put_var( ncid, dzw_varid, mesh % dzw ) )
        CALL Check( nf90_put_var( ncid, dxt_varid, mesh % dxt ) )
        CALL Check( nf90_put_var( ncid, dyt_varid, mesh % dyt ) )
        CALL Check( nf90_put_var( ncid, tarea_varid, mesh % tArea ) )
        CALL Check( nf90_put_var( ncid, tmask_varid, mesh % tracermask ) )
        PRINT*,'                            Done!' 
        ! Close the netcdf file
        CALL Check( nf90_close( ncid ) )

 END SUBROUTINE WriteNetCDF_POP_Mesh
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
      theGrid % tracerMask = 0.0_prec

      ! Count the number of wet-points and set the IJK to DOF mapping
      DO j = 1, theGrid % nY
         DO i = 1, theGrid % nX

            IF( theGrid % KMT(i,j) > 0 )THEN
              DO k = 1, theGrid % KMT(i,j)
                 nDOF = nDOF + 1
                 theGrid % tracerMask(i,j,k) = 1.0_prec
                 theGrid % IJKtoDOF(i,j,k)   = nDOF
              ENDDO
            ENDIF

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

            IF( theGrid % KMT(i,j) > 0 )THEN
              DO k = 1, theGrid % KMT(i,j)
                 nDOF = nDOF + 1
                 theGrid % DOFtoIJK(1,nDOF) = i
                 theGrid % DOFtoIJK(2,nDOF) = j
                 theGrid % DOFtoIJK(3,nDOF) = k
              ENDDO
            ENDIF

         ENDDO
      ENDDO

 END SUBROUTINE ConstructWetPointMap_POP_Mesh
!
 SUBROUTINE ConstructWetPointMapWithMask_POP_Mesh( theGrid )
 ! S/R ConstructWetPointMapWithMask
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

      ! Count the number of wet-points and set the IJK to DOF mapping
      DO j = 1, theGrid % nY
         DO i = 1, theGrid % nX
            DO k = 1, theGrid % KMT(i,j)
               IF( theGrid % tracerMask(i,j,k) /= fillValue )THEN
                  nDOF = nDOF + 1
                  theGrid % IJKtoDOF(i,j,k)   = nDOF
               ENDIF
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
               IF( theGrid % tracerMask(i,j,k) /= fillValue )THEN
                  nDOF = nDOF + 1
                  theGrid % DOFtoIJK(1,nDOF) = i
                  theGrid % DOFtoIJK(2,nDOF) = j
                  theGrid % DOFtoIJK(3,nDOF) = k
               ENDIF
            ENDDO
         ENDDO
      ENDDO

 END SUBROUTINE ConstructWetPointMapWithMask_POP_Mesh
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

         i = theGrid % DOFtoIJK(1,l)
         j = theGrid % DOFtoIJK(2,l)
         k = theGrid % DOFtoIJK(3,l)

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
   REAL(prec)        :: ijkArray(1:theGrid % nX, 1:theGrid % nY, 1:theGrid % nZ)
   REAL(prec)        :: dofArray(1:theGrid % nDOF)
   ! Local
   INTEGER :: i, j, k, l

      DO l = 1, theGrid % nDof

         i = theGrid % DOFtoIJK(1,l)
         j = theGrid % DOFtoIJK(2,l)
         k = theGrid % DOFtoIJK(3,l)

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
END MODULE POP_Mesh_Class
