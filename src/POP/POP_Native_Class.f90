! POP_Native_Class.f90
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
      REAL(prec), ALLOCATABLE  :: mask(:,:,:,:)
      REAL(prec), ALLOCATABLE  :: source(:,:,:,:)
      REAL(prec), ALLOCATABLE  :: rFac(:,:,:,:)

      REAL(prec), ALLOCATABLE  :: volume(:,:,:)
      REAL(prec), ALLOCATABLE  :: temperature(:,:,:)
      REAL(prec), ALLOCATABLE  :: salinity(:,:,:)
      REAL(prec), ALLOCATABLE  :: density(:,:,:)

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

      PROCEDURE :: LoadOceanState           => LoadOceanState_POP_Native
      PROCEDURE :: WriteOceanState          => WriteOceanState_POP_Native

   END TYPE POP_Native


   INTEGER :: ncid_PN
   INTEGER :: z_dimid_PN, y_dimid_PN, x_dimid_PN, rec_dimid_PN
   INTEGER :: z_varid_PN, y_varid_PN, x_varid_PN
   INTEGER, ALLOCATABLE :: tracer_varid_PN(:)
   INTEGER, ALLOCATABLE :: source_varid_PN(:)
   INTEGER, ALLOCATABLE :: rfac_varid_PN(:)
   INTEGER, ALLOCATABLE :: mask_varid_PN(:)
   INTEGER              :: volume_varid_PN

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
                this % mask(1:nX,1:nY,1:nZ,1:nTracers), &
                this % source(1:nX,1:nY,1:nZ,1:nTracers), &
                this % rFac(1:nX,1:nY,1:nZ,1:nTracers), &
                this % volume(1:nX,1:nY,1:nZ), &
                this % temperature(1:nX,1:nY,1:nZ), &
                this % salinity(1:nX,1:nY,1:nZ), &
                this % density(1:nX,1:nY,1:nZ) )

      this % tracer       = 0.0_prec
      this % mask         = 1.0_prec
      this % source       = 0.0_prec
      this % rFac         = 0.0_prec
      this % volume       = 0.0_prec
      this % temperature  = 0.0_prec
      this % salinity     = 0.0_prec
      this % density     = 0.0_prec

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
                  this % mask, &
                  this % source, &
                  this % rFac, &
                  this % volume, &
                  this % temperature, &
                  this % salinity, &
                  this % density )

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
                                      volume_varid_PN ) )

         CALL Check( nf90_def_var( ncid_PN, "Source_Particulate", NF90_DOUBLE, &
                                      (/ x_dimid_PN, y_dimid_PN, z_dimid_PN /), &
                                      source_varid_PN(1) ) )
         CALL Check( nf90_def_var( ncid_PN, "rFac_Particulate", NF90_DOUBLE, &
                                      (/ x_dimid_PN, y_dimid_PN, z_dimid_PN /), &
                                      rFac_varid_PN(1) ) )
         CALL Check( nf90_def_var( ncid_PN, "mask_Particulate", NF90_DOUBLE, &
                                      (/ x_dimid_PN, y_dimid_PN, z_dimid_PN /), &
                                      mask_varid_PN(1) ) )

         CALL Check( nf90_def_var( ncid_PN, "Source_Radionuclide", NF90_DOUBLE, &
                                      (/ x_dimid_PN, y_dimid_PN, z_dimid_PN /), &
                                      source_varid_PN(2) ) )
         CALL Check( nf90_def_var( ncid_PN, "rFac_Radionuclide", NF90_DOUBLE, &
                                      (/ x_dimid_PN, y_dimid_PN, z_dimid_PN /), &
                                      rFac_varid_PN(2) ) )
         CALL Check( nf90_def_var( ncid_PN, "mask_Radionuclide", NF90_DOUBLE, &
                                      (/ x_dimid_PN, y_dimid_PN, z_dimid_PN /), &
                                      mask_varid_PN(2) ) )

      ELSEIF( modelType == DyeModel .OR. modelType == SettlingModel )THEN
      
         DO i = 1, this % nTracers
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

            CALL Check( nf90_put_att( ncid_PN, tracer_varid_PN(i), "long_name", &
                                      "Dye tracer concentration of tracer "//tracerid ) )
            CALL Check( nf90_put_att( ncid_PN, tracer_varid_PN(i), "units", "" ) )
            CALL Check( nf90_put_att( ncid_PN, tracer_varid_PN(i), "coordinates", &
                                      "TLONG TLAT z_t" ) )
            CALL Check( nf90_put_att( ncid_PN, tracer_varid_PN(i), "_FillValue", &
                                      fillValue) )
            CALL Check( nf90_put_att( ncid_PN, tracer_varid_PN(i), "missing_value", &
                                      fillValue) )

            CALL Check( nf90_put_att( ncid_PN, source_varid_PN(i), "long_name", &
                                      "Relaxation field of tracer "//tracerid ) )
            CALL Check( nf90_put_att( ncid_PN, source_varid_PN(i), "units", "" ) )
            CALL Check( nf90_put_att( ncid_PN, source_varid_PN(i), "coordinates", &
                                      "TLONG TLAT z_t" ) )
            CALL Check( nf90_put_att( ncid_PN, source_varid_PN(i), "_FillValue", &
                                      fillValue) )
            CALL Check( nf90_put_att( ncid_PN, source_varid_PN(i), "missing_value", &
                                      fillValue) )

            CALL Check( nf90_put_att( ncid_PN, rFac_varid_PN(i), "long_name", &
                                      "Relaxation frequency of tracer "//tracerid ) )
            CALL Check( nf90_put_att( ncid_PN, rFac_varid_PN(i), "units", "" ) )
            CALL Check( nf90_put_att( ncid_PN, rFac_varid_PN(i), "coordinates", &
                                      "TLONG TLAT z_t" ) )
            CALL Check( nf90_put_att( ncid_PN, rFac_varid_PN(i), "_FillValue", &
                                      fillValue) )
            CALL Check( nf90_put_att( ncid_PN, rFac_varid_PN(i), "missing_value", &
                                      fillValue) )

            CALL Check( nf90_put_att( ncid_PN, mask_varid_PN(i), "long_name", &
                                      "Source-region mask of tracer "//tracerid ) )
            CALL Check( nf90_put_att( ncid_PN, mask_varid_PN(i), "units", "" ) )
            CALL Check( nf90_put_att( ncid_PN, mask_varid_PN(i), "coordinates", &
                                      "TLONG TLAT z_t" ) )
            CALL Check( nf90_put_att( ncid_PN, mask_varid_PN(i), "_FillValue", &
                                      fillValue) )
            CALL Check( nf90_put_att( ncid_PN, mask_varid_PN(i), "missing_value", &
                                      fillValue) )

         ENDDO
      
            CALL Check( nf90_def_var( ncid_PN, "VolumeCorrection", NF90_DOUBLE, &
                                      (/ x_dimid_PN, y_dimid_PN, z_dimid_PN, rec_dimid_PN /), &
                                      volume_varid_PN )  )

            CALL Check( nf90_put_att( ncid_PN, volume_varid_PN, "long_name", &
                                      "Fractional change of fluid volume" ) )
            CALL Check( nf90_put_att( ncid_PN, volume_varid_PN, "units", "" ) )
            CALL Check( nf90_put_att( ncid_PN, volume_varid_PN, "coordinates", &
                                      "TLONG TLAT z_t" ) )
            CALL Check( nf90_put_att( ncid_PN, volume_varid_PN, "_FillValue", &
                                      fillValue) )
            CALL Check( nf90_put_att( ncid_PN, volume_varid_PN, "missing_value", &
                                      fillValue) )


      ENDIF

      ! And assign attributes to each variable

      CALL Check( nf90_put_att( ncid_PN, x_varid_PN, "long_name", "array of t-grid longitudes" ) )
      CALL Check( nf90_put_att( ncid_PN, x_varid_PN, "units", "degrees_east" ) )
      CALL Check( nf90_put_att( ncid_PN, x_varid_PN, "_FillValue", fillValue) )
      CALL Check( nf90_put_att( ncid_PN, x_varid_PN, "missing_value", fillValue) )

      CALL Check( nf90_put_att( ncid_PN, y_varid_PN, "long_name", "array of t-grid latitudes" ) )
      CALL Check( nf90_put_att( ncid_PN, y_varid_PN, "units", "degrees_north" ) )
      CALL Check( nf90_put_att( ncid_PN, y_varid_PN, "_FillValue", fillValue) )
      CALL Check( nf90_put_att( ncid_PN, y_varid_PN, "missing_value", fillValue) )

      CALL Check( nf90_put_att( ncid_PN, z_varid_PN, "long_name", "depth from surface to midpoint of layer" ) )
      CALL Check( nf90_put_att( ncid_PN, z_varid_PN, "units", "centimeters" ) )
      CALL Check( nf90_put_att( ncid_PN, z_varid_PN, "positive", "down" ) )


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
                                     volume_varid_PN ) )
      
         CALL Check( nf90_inq_varid( ncid_PN, "Source_Particulate", &
                                      source_varid_PN(1) ) )
         CALL Check( nf90_inq_varid( ncid_PN, "rFac_Particulate", &
                                      rFac_varid_PN(1) ) )
         CALL Check( nf90_inq_varid( ncid_PN, "mask_Particulate", &
                                      mask_varid_PN(1) ) )
         
         CALL Check( nf90_inq_varid( ncid_PN, "Source_Radionuclide", &
                                      source_varid_PN(2) ) )
         CALL Check( nf90_inq_varid( ncid_PN, "rFac_Radionuclide", &
                                      rFac_varid_PN(2) ) )
         CALL Check( nf90_inq_varid( ncid_PN, "mask_Radionuclide", &
                                      mask_varid_PN(2) ) )

      ELSEIF( modelType == DyeModel .OR. modelType == SettlingModel )THEN
      
         DO i = 1, this % nTracers
            WRITE( tracerid, '(I2.2)') i
            CALL Check( nf90_inq_varid( ncid_PN, "DyeTracer_"//tracerid, &
                                      tracer_varid_PN(i) ) )
            CALL Check( nf90_inq_varid( ncid_PN, "Source_"//tracerid, &
                                      source_varid_PN(i) ) )
            CALL Check( nf90_inq_varid( ncid_PN, "rFac_"//tracerid, &
                                      rFac_varid_PN(i) ) )
            CALL Check( nf90_inq_varid( ncid_PN, "mask_"//tracerid, &
                                      mask_varid_PN(i) ) )
         ENDDO
            CALL Check( nf90_inq_varid( ncid_PN, "VolumeCorrection", &
                                      volume_varid_PN ) )
      
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
       
 END SUBROUTINE FinalizeNetCDF_POP_Native
!
 SUBROUTINE LoadOceanState_POP_Native( this, mesh, filename )
   IMPLICIT NONE
   CLASS( POP_Native ), INTENT(inout) :: this
   TYPE( POP_Mesh ), INTENT(in)       :: mesh
   CHARACTER(*), INTENT(in)           :: filename
   ! Local
   INTEGER :: ncid, varid
   INTEGER :: start2D(1:2), recCount2D(1:2)
   INTEGER :: start3D(1:3), recCount3D(1:3)

      CALL Check( nf90_open( TRIM(filename), nf90_nowrite, ncid ) )

      start2D    = (/1, 1/)
      recCount2D = (/mesh % nX, mesh % nY/)
      start3D    = (/1, 1, 1/)
      recCount3D = (/mesh % nX, mesh % nY, mesh % nZ/)

      PRINT*, 'Loading SSH'
      CALL Check( nf90_inq_varid( ncid, "SSH",varid ) )
      CALL Check( nf90_get_var( ncid, &
                                varid, &
                                this % volume(:,:,1), &
                                start2D, recCount2D ) )

      PRINT*, 'Loading TEMP'
      CALL Check( nf90_inq_varid( ncid, "TEMP",varid ) )
      CALL Check( nf90_get_var( ncid, &
                                varid, &
                                this % temperature, &
                                start3D, recCount3D ) )

      PRINT*, 'Loading SALT'
      CALL Check( nf90_inq_varid( ncid, "SALT",varid ) )
      CALL Check( nf90_get_var( ncid, &
                                varid, &
                                this % salinity, &
                                start3D, recCount3D ) )

      PRINT*, 'Loading PD'
      CALL Check( nf90_inq_varid( ncid, "PD",varid ) )
      CALL Check( nf90_get_var( ncid, &
                                varid, &
                                this % density, &
                                start3D, recCount3D ) )

      PRINT*, 'DONE'

      CALL Check( nf90_close( ncid ) )

 END SUBROUTINE LoadOceanState_POP_Native
!
 SUBROUTINE WriteOceanState_POP_Native( this, mesh, filename )
   IMPLICIT NONE
   CLASS( POP_Native ), INTENT(inout) :: this
   TYPE( POP_Mesh ), INTENT(in)       :: mesh
   CHARACTER(*), INTENT(in)           :: filename
   ! Local
   INTEGER :: ncid, z_dimid, x_dimid, y_dimid
   INTEGER :: varid_ssh, varid_temp, varid_salt, varid_pd
   INTEGER :: start2D(1:2), recCount2D(1:2)
   INTEGER :: start3D(1:3), recCount3D(1:3)

      CALL Check( nf90_create( PATH=TRIM(filename),&
                               CMODE=OR(nf90_clobber,nf90_64bit_offset),&
                               NCID=ncid ) )
      ! Create the dimensions - the dimension names are currently chosen based
      CALL Check( nf90_def_dim( ncid, "z_t", mesh % nZ, z_dimid ) ) 
      CALL Check( nf90_def_dim( ncid, "nlon", mesh % nX, x_dimid ) ) 
      CALL Check( nf90_def_dim( ncid, "nlat", mesh % nY, y_dimid ) ) 


      start2D    = (/1, 1/)
      recCount2D = (/mesh % nX, mesh % nY/)
      start3D    = (/1, 1, 1/)
      recCount3D = (/mesh % nX, mesh % nY, mesh % nZ/)

      CALL Check( nf90_def_var( ncid, "SSH", NF90_DOUBLE,&
                                (/ x_dimid, y_dimid /), &
                                 varid_ssh ) )

      CALL Check( nf90_def_var( ncid, "TEMP", NF90_DOUBLE,&
                                (/ x_dimid, y_dimid, z_dimid /), &
                                 varid_temp ) )

      CALL Check( nf90_def_var( ncid, "SALT", NF90_DOUBLE,&
                                (/ x_dimid, y_dimid, z_dimid /), &
                                 varid_salt ) )

      CALL Check( nf90_def_var( ncid, "PD", NF90_DOUBLE,&
                                (/x_dimid, y_dimid, z_dimid/),&
                                 varid_pd ) )


      CALL Check( nf90_enddef(ncid) )

      CALL Check( nf90_put_var( ncid, &
                                varid_ssh, &
                                this % volume(:,:,1), &
                                start2D, recCount2D ) )

      CALL Check( nf90_put_var( ncid, &
                                varid_temp, &
                                this % temperature, &
                                start3D, recCount3D ) )

      CALL Check( nf90_put_var( ncid, &
                                varid_salt, &
                                this % salinity, &
                                start3D, recCount3D ) )

      CALL Check( nf90_put_var( ncid, &
                                varid_pd, &
                                this % density, &
                                start3D, recCount3D ) )
      CALL Check( nf90_close( ncid ) )

 END SUBROUTINE WriteOceanState_POP_Native
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

      DO i = 1, this % nTracers
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
         CALL Check( nf90_put_var( ncid_PN, &
                                   volume_varid_pn, &
                                   this % volume, &
                                   start, recCount ) )      

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

         DO i = 1, this % nTracers
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
            CALL Check( nf90_get_var( ncid_PN, &
                                      volume_varid_pn, &
                                      this % volume, &
                                      start, recCount ) )

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
