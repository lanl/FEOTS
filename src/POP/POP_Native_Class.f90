! POP_Native_Class.f90
! 
! Copyright 2016 Joseph Schoonover, Los Alamos National Laboratory (jschoonover@lanl.gov) 
! All rights reserved. 
! 
! POP_Native_Class.f90 is part of the Fast Equilibration of Ocean Tracers Software (FEOTS). 
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !
 
 
 MODULE POP_Native_Class

! src/common/
 USE ModelPrecision
 USE CommonRoutines
! src/POP/
 USE POP_Mesh_Class
!
USE netcdf


 IMPLICIT NONE


   TYPE POP_Native
      INTEGER                  :: nX, nY, nZ, nTracers
      REAL(prec), ALLOCATABLE  :: tracer(:,:,:,:)
      REAL(prec), ALLOCATABLE  :: hardSet(:,:,:,:)
      REAL(prec), ALLOCATABLE  :: mask(:,:,:,:)
      REAL(prec), ALLOCATABLE  :: source(:,:,:,:)
      REAL(prec), ALLOCATABLE  :: rFac(:,:,:,:)

      CONTAINS

      PROCEDURE :: Build => Build_POP_Native
      PROCEDURE :: Trash => Trash_POP_Native

      PROCEDURE :: InitializeForNetCDFRead  => InitializeForNetCDFRead_POP_Native
      PROCEDURE :: InitializeForNetCDFWrite => InitializeForNetCDFWrite_POP_Native
      PROCEDURE :: FinalizeNetCDF           => FinalizeNetCDF_POP_Native
      PROCEDURE :: WriteNetCDFRecord        => WriteNetCDFRecord_POP_Native
      PROCEDURE :: ReadNetCDFRecord         => ReadNetCDFRecord_POP_Native
      PROCEDURE :: WriteSourceEtcNetCDF     => WriteSourceEtcNetCDF_POP_Native
      PROCEDURE :: ReadSourceEtcNetCDF      => ReadSourceEtcNetCDF_POP_Native
      PROCEDURE :: LoadTracerFromNetCDF     => LoadTracerFromNetCDF_POP_Native

   END TYPE POP_Native


   INTEGER :: ncid_PN
   INTEGER :: z_dimid_PN, y_dimid_PN, x_dimid_PN, rec_dimid_PN
   INTEGER :: z_varid_PN, y_varid_PN, x_varid_PN
   INTEGER, ALLOCATABLE :: tracer_varid_PN(:)
   INTEGER, ALLOCATABLE :: source_varid_PN(:)
   INTEGER, ALLOCATABLE :: rfac_varid_PN(:)
   INTEGER, ALLOCATABLE :: mask_varid_PN(:)
   INTEGER, ALLOCATABLE :: hardset_varid_PN(:)
CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE Build_POP_Native( this, mesh, nTracers ) 
 ! S/R Build
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( POP_Native ), INTENT(out) :: this
   TYPE( POP_Mesh ), INTENT(in)     :: mesh
   INTEGER, INTENT(in)              :: nTracers
   ! Local 
   INTEGER :: nX, nY, nZ

      PRINT*, 'S/R : Build_POP_Native : Start...'
      nX = mesh % nX
      nY = mesh % nY
      nZ = mesh % nZ
      this % nTracers = nTracers
      this % nX = nX
      this % nY = nY
      this % nZ = nZ

      ALLOCATE( this % tracer(1:nX,1:nY,1:nZ,1:nTracers), &
                this % hardSet(1:nX,1:nY,1:nZ,1:nTracers), &
                this % mask(1:nX,1:nY,1:nZ,1:nTracers), &
                this % source(1:nX,1:nY,1:nZ,1:nTracers), &
                this % rFac(1:nX,1:nY,1:nZ,1:nTracers) )
      this % tracer  = 0.0_prec
      this % hardSet = 0.0_prec
      this % mask    = 1.0_prec
      this % source  = 0.0_prec
      this % rFac    = 0.0_prec

      PRINT*, 'S/R : Build_POP_Native : Finish.'

 END SUBROUTINE Build_POP_Native
!
 SUBROUTINE Trash_POP_Native( this )
 ! S/R Trash
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( POP_Native ), INTENT(inout) :: this

      DEALLOCATE( this % tracer, &
                  this % hardSet, &
                  this % mask, &
                  this % source, &
                  this % rFac )

 END SUBROUTINE Trash_POP_Native
!
 SUBROUTINE InitializeForNetCDFWrite_POP_Native( this, modelType, mesh, filename )
   IMPLICIT NONE
   CLASS( POP_Native ), INTENT(in) :: this
   INTEGER, INTENT(in)             :: modelType
   TYPE( POP_Mesh ),  INTENT(in)   :: mesh
   CHARACTER(*), INTENT(in)        :: filename
   ! Local
   INTEGER :: i
   CHARACTER(2) :: tracerid

      ! Create the netcdf file and generate a file handle referenced by the
      ! integer "ncid_PN"

      CALL Check( nf90_create( PATH=TRIM(filename),&
                               CMODE=OR(nf90_clobber,nf90_64bit_offset),&
                               NCID=ncid_PN ) )
      ! Create the dimensions - the dimension names are currently chosen based
      ! on the netcdf output from POP (version ???, branch ???)
      CALL Check( nf90_def_dim( ncid_PN, "z_t", this % nZ, z_dimid_PN ) ) 
      CALL Check( nf90_def_dim( ncid_PN, "nlon", this % nX, x_dimid_PN ) ) 
      CALL Check( nf90_def_dim( ncid_PN, "nlat", this % nY, y_dimid_PN ) ) 
      CALL Check( nf90_def_dim( ncid_PN, "time", NF90_UNLIMITED, rec_dimid_PN ) )

      ! Create variables -- here we need to create arrays for the dimensions
      CALL Check( nf90_def_var( ncid_PN, "z_t", NF90_DOUBLE, z_dimid_PN, z_varid_PN ) )
      CALL Check( nf90_def_var( ncid_PN, "TLAT", NF90_DOUBLE, (/ x_dimid_PN, y_dimid_PN /), y_varid_PN ) )
      CALL Check( nf90_def_var( ncid_PN, "TLONG", NF90_DOUBLE, (/ x_dimid_PN, y_dimid_PN /), x_varid_PN ) )

      ALLOCATE( tracer_varid_PN(1:this % nTracers) )
      ALLOCATE( source_varid_PN(1:this % nTracers) )
      ALLOCATE( rfac_varid_PN(1:this % nTracers) )
      ALLOCATE( mask_varid_PN(1:this % nTracers) )
      ALLOCATE( hardset_varid_PN(1:this % nTracers) )
      ! Set up the tracer field names based on the model type
      IF( modelType == ImpulseField )THEN
      
         ! Need to check with Wilbert if this is the correct naming 
         ! convention for reporting the impulse fields.
         DO i = 1, this % nTracers
            WRITE( tracerid, '(I2.2)') i
            CALL Check( nf90_def_var( ncid_PN, "ADV_3D_IRF_"//tracerid, NF90_DOUBLE,&
                                      (/ x_dimid_PN, y_dimid_PN, z_dimid_PN /), &
                                      tracer_varid_PN(i) ) )
         ENDDO
         
      
      ELSEIF( modelType == RadionuclideModel )THEN
      
         CALL Check( nf90_def_var( ncid_PN, "Particulate", NF90_DOUBLE, &
                                      (/ x_dimid_PN, y_dimid_PN, z_dimid_PN, rec_dimid_PN /), &
                                      tracer_varid_PN(1) ) )
                                      
         CALL Check( nf90_def_var( ncid_PN, "Radionuclide", NF90_DOUBLE, &
                                      (/ x_dimid_PN, y_dimid_PN, z_dimid_PN, rec_dimid_PN /), &
                                      tracer_varid_PN(2) ) )

         CALL Check( nf90_def_var( ncid_PN, "VolumeCorrection", NF90_DOUBLE, &
                                      (/ x_dimid_PN, y_dimid_PN, z_dimid_PN, rec_dimid_PN /), &
                                      tracer_varid_PN(3) ) )

         CALL Check( nf90_def_var( ncid_PN, "Source_Particulate", NF90_DOUBLE, &
                                      (/ x_dimid_PN, y_dimid_PN, z_dimid_PN /), &
                                      source_varid_PN(1) ) )
         CALL Check( nf90_def_var( ncid_PN, "rFac_Particulate", NF90_DOUBLE, &
                                      (/ x_dimid_PN, y_dimid_PN, z_dimid_PN /), &
                                      rFac_varid_PN(1) ) )
         CALL Check( nf90_def_var( ncid_PN, "mask_Particulate", NF90_DOUBLE, &
                                      (/ x_dimid_PN, y_dimid_PN, z_dimid_PN /), &
                                      mask_varid_PN(1) ) )
         CALL Check( nf90_def_var( ncid_PN, "hardset_Particulate", NF90_DOUBLE, &
                                      (/ x_dimid_PN, y_dimid_PN, z_dimid_PN /), &
                                      hardset_varid_PN(1) ) )


         CALL Check( nf90_def_var( ncid_PN, "Source_Radionuclide", NF90_DOUBLE, &
                                      (/ x_dimid_PN, y_dimid_PN, z_dimid_PN /), &
                                      source_varid_PN(2) ) )
         CALL Check( nf90_def_var( ncid_PN, "rFac_Radionuclide", NF90_DOUBLE, &
                                      (/ x_dimid_PN, y_dimid_PN, z_dimid_PN /), &
                                      rFac_varid_PN(2) ) )
         CALL Check( nf90_def_var( ncid_PN, "mask_Radionuclide", NF90_DOUBLE, &
                                      (/ x_dimid_PN, y_dimid_PN, z_dimid_PN /), &
                                      mask_varid_PN(2) ) )
         CALL Check( nf90_def_var( ncid_PN, "hardset_Radionuclide", NF90_DOUBLE, &
                                      (/ x_dimid_PN, y_dimid_PN, z_dimid_PN /), &
                                      hardset_varid_PN(2) ) )

      ELSEIF( modelType == DyeModel .OR. modelType == SettlingModel )THEN
      
         DO i = 1, this % nTracers-1
            WRITE( tracerid, '(I2.2)') i
            CALL Check( nf90_def_var( ncid_PN, "DyeTracer_"//tracerid, NF90_DOUBLE, &
                                      (/ x_dimid_PN, y_dimid_PN, z_dimid_PN, rec_dimid_PN /), &
                                      tracer_varid_PN(i) ) )
            CALL Check( nf90_def_var( ncid_PN, "Source_"//tracerid, NF90_DOUBLE, &
                                      (/ x_dimid_PN, y_dimid_PN, z_dimid_PN /), &
                                      source_varid_PN(i) ) )
            CALL Check( nf90_def_var( ncid_PN, "rFac_"//tracerid, NF90_DOUBLE, &
                                      (/ x_dimid_PN, y_dimid_PN, z_dimid_PN /), &
                                      rFac_varid_PN(i) ) )
            CALL Check( nf90_def_var( ncid_PN, "mask_"//tracerid, NF90_DOUBLE, &
                                      (/ x_dimid_PN, y_dimid_PN, z_dimid_PN /), &
                                      mask_varid_PN(i) ) )
            CALL Check( nf90_def_var( ncid_PN, "hardset_"//tracerid, NF90_DOUBLE, &
                                      (/ x_dimid_PN, y_dimid_PN, z_dimid_PN /), &
                                      hardset_varid_PN(i) ) )
         ENDDO
      
            CALL Check( nf90_def_var( ncid_PN, "VolumeCorrection", NF90_DOUBLE, &
                                      (/ x_dimid_PN, y_dimid_PN, z_dimid_PN, rec_dimid_PN /), &
                                      tracer_varid_PN(this%nTracers) ) )

      ENDIF

      ! And assign units to each variable
      CALL Check( nf90_put_att( ncid_PN, z_varid_PN, "units", "cm" ) )
      CALL Check( nf90_put_att( ncid_PN, y_varid_PN, "units", "deg N" ) )
      CALL Check( nf90_put_att( ncid_PN, z_varid_PN, "units", "deg E" ) )

      ! End the Define Mode
      CALL Check( nf90_enddef(ncid_PN) )
      ! 
      CALL Check( nf90_put_var( ncid_PN, z_varid_PN, mesh % z ) )
      CALL Check( nf90_put_var( ncid_PN, y_varid_PN, mesh % tLat ) )
      CALL Check( nf90_put_var( ncid_PN, x_varid_PN, mesh % tLon ) )


 END SUBROUTINE InitializeForNetCDFWrite_POP_Native
! 
 SUBROUTINE InitializeForNetCDFRead_POP_Native( this, modelType, filename )
   IMPLICIT NONE
   CLASS( POP_Native ), INTENT(in) :: this
   INTEGER, INTENT(in)             :: modelType
   CHARACTER(*), INTENT(in)        :: filename
   ! Local
   INTEGER :: i
   CHARACTER(2) :: tracerid

      ! Create the netcdf file and generate a file handle referenced by the
      ! integer "ncid_PN"
      CALL Check( nf90_open( TRIM(filename), nf90_nowrite, ncid_PN ) )
      ! Create variables -- here we need to create arrays for the dimensions
      CALL Check( nf90_inq_varid( ncid_PN, "z_t", z_varid_PN ) )
      CALL Check( nf90_inq_varid( ncid_PN, "TLAT", y_varid_PN ) )
      CALL Check( nf90_inq_varid( ncid_PN, "TLONG", x_varid_PN ) )

      ALLOCATE( tracer_varid_PN(1:this % nTracers) )
      ALLOCATE( source_varid_PN(1:this % nTracers) )
      ALLOCATE( rfac_varid_PN(1:this % nTracers) )
      ALLOCATE( mask_varid_PN(1:this % nTracers) )
      ALLOCATE( hardset_varid_PN(1:this % nTracers) )
      ! Set up the tracer field names based on the model type
      IF( modelType == ImpulseResponseField )THEN
      
         ! Need to check with Wilbert if this is the correct naming 
         ! convention for reporting the impulse fields.
         DO i = 1, this % nTracers-1
            WRITE( tracerid, '(I2.2)') i
            PRINT*, "ADV_3D_IRF_"//tracerid
            CALL Check( nf90_inq_varid( ncid_PN, "ADV_3D_IRF_"//tracerid, &
                                        tracer_varid_PN(i) ) )
         ENDDO
         CALL Check( nf90_inq_varid( ncid_PN, "VDC_S", &
                                     tracer_varid_PN(this % nTracers) ) )
         PRINT*, "VDC_S"

      ELSEIF( modelType == ImpulseField )THEN
      
         DO i = 1, this % nTracers
            WRITE( tracerid, '(I2.2)') i
            CALL Check( nf90_inq_varid( ncid_PN, "ADV_3D_IRF_"//tracerid, &
                                        tracer_varid_PN(i) ) )
         ENDDO
                             
      
      ELSEIF( modelType == RadionuclideModel )THEN
      
         CALL Check( nf90_inq_varid( ncid_PN, "Particulate", &
                                      tracer_varid_PN(1) ) )
                                      
         CALL Check( nf90_inq_varid( ncid_PN, "Radionuclide", &
                                      tracer_varid_PN(2) ) )
      
         CALL Check( nf90_inq_varid( ncid_PN, "VolumeCorrection", &
                                      tracer_varid_PN(3) ) )
      
         CALL Check( nf90_inq_varid( ncid_PN, "Source_Particulate", &
                                      source_varid_PN(1) ) )
         CALL Check( nf90_inq_varid( ncid_PN, "rFac_Particulate", &
                                      rFac_varid_PN(1) ) )
         CALL Check( nf90_inq_varid( ncid_PN, "mask_Particulate", &
                                      mask_varid_PN(1) ) )
         CALL Check( nf90_inq_varid( ncid_PN, "hardset_Particulate", &
                                      hardset_varid_PN(1) ) )
         
         CALL Check( nf90_inq_varid( ncid_PN, "Source_Radionuclide", &
                                      source_varid_PN(2) ) )
         CALL Check( nf90_inq_varid( ncid_PN, "rFac_Radionuclide", &
                                      rFac_varid_PN(2) ) )
         CALL Check( nf90_inq_varid( ncid_PN, "mask_Radionuclide", &
                                      mask_varid_PN(2) ) )
         CALL Check( nf90_inq_varid( ncid_PN, "hardset_Radionuclide", &
                                      hardset_varid_PN(2) ) )

      ELSEIF( modelType == DyeModel .OR. modelType == SettlingModel )THEN
      
         DO i = 1, this % nTracers-1
            WRITE( tracerid, '(I2.2)') i
            CALL Check( nf90_inq_varid( ncid_PN, "DyeTracer_"//tracerid, &
                                      tracer_varid_PN(i) ) )
            CALL Check( nf90_inq_varid( ncid_PN, "Source_"//tracerid, &
                                      source_varid_PN(i) ) )
            CALL Check( nf90_inq_varid( ncid_PN, "rFac_"//tracerid, &
                                      rFac_varid_PN(i) ) )
            CALL Check( nf90_inq_varid( ncid_PN, "mask_"//tracerid, &
                                      mask_varid_PN(i) ) )
            CALL Check( nf90_inq_varid( ncid_PN, "hardset_"//tracerid, &
                                      hardset_varid_PN(i) ) )
         ENDDO
            CALL Check( nf90_inq_varid( ncid_PN, "VolumeCorrection", &
                                      tracer_varid_PN(this % nTracers) ) )
      
      ENDIF

 END SUBROUTINE InitializeForNetCDFRead_POP_Native
!   
 SUBROUTINE FinalizeNetCDF_POP_Native( this )
    IMPLICIT NONE
    CLASS( POP_Native ) :: this
 
       CALL Check( nf90_close( ncid_PN ) )
       DEALLOCATE( tracer_varid_PN )
       DEALLOCATE( source_varid_PN )
       DEALLOCATE( rFac_varid_PN )
       DEALLOCATE( mask_varid_PN )
       DEALLOCATE( hardset_varid_PN )
       
 END SUBROUTINE FinalizeNetCDF_POP_Native
!
 SUBROUTINE WriteSourceEtcNETCDF_POP_Native( this, mesh )
   IMPLICIT NONE
   CLASS( POP_Native ), INTENT(in) :: this
   TYPE( POP_Mesh ), INTENT(in)    :: mesh
   ! Local
   INTEGER :: start(1:3), recCount(1:3)
   INTEGER :: i

      start = (/1, 1, 1/)
      recCount = (/mesh % nX, mesh % nY, mesh % nZ/)

      DO i = 1, this % nTracers-1
         CALL Check( nf90_put_var( ncid_PN, &
                                   source_varid_pn(i), &
                                   this % source(:,:,:,i), &
                                   start, recCount ) )      
         CALL Check( nf90_put_var( ncid_PN, &
                                   rfac_varid_pn(i), &
                                   this % rfac(:,:,:,i), &
                                   start, recCount ) )      
         CALL Check( nf90_put_var( ncid_PN, &
                                   mask_varid_pn(i), &
                                   this % mask(:,:,:,i), &
                                   start, recCount ) )      
         CALL Check( nf90_put_var( ncid_PN, &
                                   hardset_varid_pn(i), &
                                   this % hardset(:,:,:,i), &
                                   start, recCount ) )      
      ENDDO

 END SUBROUTINE WriteSourceEtcNETCDF_POP_Native
!
 SUBROUTINE WriteNETCDFRecord_POP_Native( this, mesh, recordID )
   IMPLICIT NONE
   CLASS( POP_Native ), INTENT(in) :: this
   TYPE( POP_Mesh ), INTENT(in)    :: mesh
   INTEGER, INTENT(in)             :: recordID
   ! Local
   INTEGER :: start(1:4), recCount(1:4)
   INTEGER :: i

      start = (/1, 1, 1, recordID/)
      recCount = (/mesh % nX, mesh % nY, mesh % nZ, 1/)

      DO i = 1, this % nTracers
         CALL Check( nf90_put_var( ncid_PN, &
                                   tracer_varid_pn(i), &
                                   this % tracer(:,:,:,i), &
                                   start, recCount ) )      
      ENDDO

 END SUBROUTINE WriteNETCDFRecord_POP_Native
!
 SUBROUTINE ReadSourceEtcNETCDF_POP_Native( this, mesh )
   IMPLICIT NONE
   CLASS( POP_Native ), INTENT(inout) :: this
   TYPE( POP_Mesh ), INTENT(in)       :: mesh
   ! Local
   INTEGER :: start(1:3), recCount(1:3)
   INTEGER :: i

         start    = (/1, 1, 1/)
         recCount = (/mesh % nX, mesh % nY, mesh % nZ/)

         DO i = 1, this % nTracers-1
            CALL Check( nf90_get_var( ncid_PN, &
                                      source_varid_pn(i), &
                                      this % source(:,:,:,i), &
                                      start, recCount ) )
            CALL Check( nf90_get_var( ncid_PN, &
                                      rfac_varid_pn(i), &
                                      this % rfac(:,:,:,i), &
                                      start, recCount ) )
            CALL Check( nf90_get_var( ncid_PN, &
                                      mask_varid_pn(i), &
                                      this % mask(:,:,:,i), &
                                      start, recCount ) )
            CALL Check( nf90_get_var( ncid_PN, &
                                      hardset_varid_pn(i), &
                                      this % hardset(:,:,:,i), &
                                      start, recCount ) )
         ENDDO

 END SUBROUTINE ReadSourceEtcNETCDF_POP_Native
!
 SUBROUTINE ReadNETCDFRecord_POP_Native( this, mesh, recordID )
   IMPLICIT NONE
   CLASS( POP_Native ), INTENT(inout) :: this
   TYPE( POP_Mesh ), INTENT(in)       :: mesh
   INTEGER, INTENT(in)                :: recordID
   ! Local
   INTEGER :: start(1:4), recCount(1:4)
   INTEGER :: i

         start    = (/1, 1, 1, recordID/)
         recCount = (/mesh % nX, mesh % nY, mesh % nZ, 1/)

         DO i = 1, this % nTracers
            CALL Check( nf90_get_var( ncid_PN, &
                                      tracer_varid_pn(i), &
                                      this % tracer(:,:,:,i), &
                                      start, recCount ) )
         ENDDO

 END SUBROUTINE ReadNETCDFRecord_POP_Native
!
 SUBROUTINE LoadTracerFromNetCDF_POP_Native( this, mesh )
   IMPLICIT NONE
   CLASS( POP_Native ), INTENT(inout) :: this
   TYPE( POP_Mesh ), INTENT(in)       :: mesh
   ! Local
   INTEGER :: start(1:3), recCount(1:3)
   INTEGER :: i

         start    = (/1, 1, 1/)
         recCount = (/mesh % nX, mesh % nY, mesh % nZ/)

         DO i = 1, this % nTracers
            CALL Check( nf90_get_var( ncid_PN, &
                                      tracer_varid_pn(i), &
                                      this % tracer(:,:,:,i), &
                                      start, recCount ) )
         ENDDO

 END SUBROUTINE LoadTracerFromNetCDF_POP_Native

END MODULE POP_Native_Class
