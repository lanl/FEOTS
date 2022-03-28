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
#include "FEOTS_Macros.h"

   TYPE NativeIO
     INTEGER :: ncid
     INTEGER :: zDimID, yDimID, xDimID, recDimID
     INTEGER :: zVarID, yVarID, xVarID
     INTEGER, ALLOCATABLE :: tracerVarID(:)
     INTEGER, ALLOCATABLE :: sourceVarID(:)
     INTEGER, ALLOCATABLE :: rfacVarID(:)
     INTEGER, ALLOCATABLE :: maskVarID(:)
     INTEGER, ALLOCATABLE :: volumeVarID(:)
   END TYPE NativeIO

   TYPE POP_IRF
      INTEGER                  :: nX, nY, nZ, nIRF
      REAL(prec), ALLOCATABLE  :: irf(:,:,:,:)

      CONTAINS

      PROCEDURE :: Build => Build_POP_IRF
      PROCEDURE :: Trash => Trash_POP_IRF

      PROCEDURE :: LoadNetCDF => LoadNetCDF_POP_IRF

   END TYPE POP_IRF

   TYPE POP_Native
      INTEGER                  :: nX, nY, nZ, nTracers
      INTEGER, ALLOCATABLE     :: tracerIds(:)
      REAL(prec), ALLOCATABLE  :: tracer(:,:,:,:)
      REAL(prec), ALLOCATABLE  :: mask(:,:,:,:)
      REAL(prec), ALLOCATABLE  :: source(:,:,:,:)
      REAL(prec), ALLOCATABLE  :: rFac(:,:,:,:)

      REAL(prec), ALLOCATABLE  :: volume(:,:,:,:)
      REAL(prec), ALLOCATABLE  :: temperature(:,:,:)
      REAL(prec), ALLOCATABLE  :: salinity(:,:,:)
      REAL(prec), ALLOCATABLE  :: density(:,:,:)

      TYPE(NativeIO) :: ioVars

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
      PROCEDURE :: WriteTracerToNetCDF      => WriteTracerToNetCDF_POP_Native

      PROCEDURE :: LoadOceanState           => LoadOceanState_POP_Native
      PROCEDURE :: WriteOceanState          => WriteOceanState_POP_Native

   END TYPE POP_Native


CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE Build_POP_IRF( this, nX, nY, nZ, nIRF ) 
 ! S/R Build
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( POP_IRF ), INTENT(out) :: this
   INTEGER, INTENT(in) :: nX, nY, nZ, nIRF

      this % nX = nX
      this % nY = nY
      this % nZ = nZ
      this % nIRF = nIRF

      ! Allocate nIRF+1 : nIRF for the IRF and +1 for the vertical diffusion
      ! coefficients
      ALLOCATE( this % irf(1:nX,1:nY,1:nZ,1:nIRF+1) )

      this % irf = 0.0_prec


 END SUBROUTINE Build_POP_IRF
!
 SUBROUTINE Trash_POP_IRF( this )
 ! S/R Trash
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( POP_IRF ), INTENT(inout) :: this

      IF(ALLOCATED(this % irf)) DEALLOCATE(this % irf)

 END SUBROUTINE Trash_POP_IRF
!
 SUBROUTINE Build_POP_Native( this, mesh, nTracers, myRank, nProcs ) 
 ! S/R Build
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( POP_Native ), INTENT(out) :: this
   TYPE( POP_Mesh ), INTENT(in)     :: mesh
   INTEGER, INTENT(in)              :: nTracers
   INTEGER, INTENT(in)              :: myRank
   INTEGER, INTENT(in)              :: nProcs
   ! Local 
   INTEGER :: nX, nY, nZ, nT, remainder, i

      nX = mesh % nX
      nY = mesh % nY
      nZ = mesh % nZ
      this % nX = nX
      this % nY = nY
      this % nZ = nZ


      ! myRank varies from 0 to nProcs-1
      ! We need to share nTracers as evenly as possible across
      ! mpi ranks
      nT = nTracers/nProcs
      remainder = nTracers - nT*nProcs
 
      ! The lower bound for the tracer IDs this rank is responsible for
      IF( myRank == nProcs-1 )THEN
        this % nTracers = nT+remainder
      ELSE
        this % nTracers = nT
      ENDIF

      ALLOCATE( this % tracer(1:nX,1:nY,1:nZ,1:this % nTracers), &
                this % mask(1:nX,1:nY,1:nZ,1:this % nTracers), &
                this % source(1:nX,1:nY,1:nZ,1:this % nTracers), &
                this % rFac(1:nX,1:nY,1:nZ,1:this % nTracers), &
                this % volume(1:nX,1:nY,1:nZ,1:this % nTracers), &
                this % tracerIds(1:this % nTracers), &
                this % ioVars % tracerVarID(1:this % nTracers), &
                this % ioVars % volumeVarID(1:this % nTracers), &
                this % ioVars % sourceVarID(1:this % nTracers), &
                this % ioVars % rfacVarID(1:this % nTracers), &
                this % ioVars % maskVarID(1:this % nTracers) )

!                this % temperature(1:nX,1:nY,1:nZ), &
!                this % salinity(1:nX,1:nY,1:nZ), &
!                this % density(1:nX,1:nY,1:nZ), &

      DO i = 1, this % nTracers
        this % tracerIds(i) =  myRank*nT+i
      ENDDO

      this % tracer       = 0.0_prec
      this % mask         = 1.0_prec
      this % source       = 0.0_prec
      this % rFac         = 0.0_prec
      this % volume       = 0.0_prec
!      this % temperature  = 0.0_prec
!      this % salinity     = 0.0_prec
!      this % density     = 0.0_prec



 END SUBROUTINE Build_POP_Native
!
 SUBROUTINE Trash_POP_Native( this )
 ! S/R Trash
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( POP_Native ), INTENT(inout) :: this

      IF(ALLOCATED(this % tracer)) DEALLOCATE(this % tracer)
      IF(ALLOCATED(this % mask)) DEALLOCATE(this % mask)
      IF(ALLOCATED(this % source)) DEALLOCATE(this % source)
      IF(ALLOCATED(this % rFac)) DEALLOCATE(this % rFac)
      IF(ALLOCATED(this % volume)) DEALLOCATE(this % volume)
!      IF(ALLOCATED(this % temperature)) DEALLOCATE(this % temperature)
!      IF(ALLOCATED(this % salinity)) DEALLOCATE(this % salinity)
!      IF(ALLOCATED(this % density )) DEALLOCATE(this % density)
      IF(ALLOCATED(this % tracerIds )) DEALLOCATE(this % tracerIds)
      IF(ALLOCATED(this % ioVars % tracerVarID)) DEALLOCATE( this % ioVars % tracerVarID)
      IF(ALLOCATED(this % ioVars % volumeVarID)) DEALLOCATE( this % ioVars % volumeVarID)
      IF(ALLOCATED(this % ioVars % sourceVarID)) DEALLOCATE( this % ioVars % sourceVarID)
      IF(ALLOCATED(this % ioVars % rfacVarID)) DEALLOCATE( this % ioVars % rfacVarID)
      IF(ALLOCATED(this % ioVars % maskVarID)) DEALLOCATE( this % ioVars % maskVarID)

 END SUBROUTINE Trash_POP_Native
!
 SUBROUTINE InitializeForNetCDFWrite_POP_Native( this, modelType, mesh, filename, initOn )
#undef __FUNC__
#define __FUNC__ "InitializeForNetCDFWrite_POP_Native"
   IMPLICIT NONE
   CLASS( POP_Native ), INTENT(inout) :: this
   INTEGER, INTENT(in)             :: modelType
   TYPE( POP_Mesh ),  INTENT(in)   :: mesh
   CHARACTER(*), INTENT(in)        :: filename
   LOGICAL, INTENT(in)             :: initOn
   ! Local
   INTEGER :: i
   CHARACTER(2) :: tracerid

      ! Create the netcdf file and generate a file handle referenced by the
      INFO('Opening NetCDF file for write.  '//TRIM(filename) )
      CALL Check( nf90_create( PATH=TRIM(filename),&
                               CMODE=OR(nf90_clobber,nf90_64bit_offset),&
                               NCID=this % ioVars % ncid ) )
      ! Create the dimensions - the dimension names are currently chosen based
      ! on the netcdf output from POP (version ???, branch ???)
      CALL Check( nf90_def_dim( this % ioVars % ncid, "z_t", this % nZ, this % ioVars % zDimID ) ) 
      CALL Check( nf90_def_dim( this % ioVars % ncid, "nlon", this % nX, this % ioVars % xDimID ) ) 
      CALL Check( nf90_def_dim( this % ioVars % ncid, "nlat", this % nY, this % ioVars % yDimID ) ) 
      CALL Check( nf90_def_dim( this % ioVars % ncid, "time", NF90_UNLIMITED, this % ioVars % recDimID ) )

      ! Create variables -- here we need to create arrays for the dimensions
      CALL Check( nf90_def_var( this % ioVars % ncid, "z_t", NF90_FLOAT, this % ioVars % zDimID, this % ioVars % zVarID ) )
      CALL Check( nf90_put_att( this % ioVars % ncid, this % ioVars % zVarID, "long_name", "depth from surface to midpoint of layer" ) )
      CALL Check( nf90_put_att( this % ioVars % ncid, this % ioVars % zVarID, "units", "centimeters" ) )
      CALL Check( nf90_put_att( this % ioVars % ncid, this % ioVars % zVarID, "positive", "down" ) )

      CALL Check( nf90_def_var( this % ioVars % ncid, "TLAT", NF90_FLOAT, (/ this % ioVars % xDimID, this % ioVars % yDimID /), this % ioVars % yVarID ) )
      CALL Check( nf90_put_att( this % ioVars % ncid, this % ioVars % yVarID, "long_name", "array of t-grid latitudes" ) )
      CALL Check( nf90_put_att( this % ioVars % ncid, this % ioVars % yVarID, "units", "degrees_north" ) )
      CALL Check( nf90_put_att( this % ioVars % ncid, this % ioVars % yVarID, "_FillValue", fillValue) )
      CALL Check( nf90_put_att( this % ioVars % ncid, this % ioVars % yVarID, "missing_value", fillValue) )

      CALL Check( nf90_def_var( this % ioVars % ncid, "TLONG", NF90_FLOAT, (/ this % ioVars % xDimID, this % ioVars % yDimID /), this % ioVars % xVarID ) )
      CALL Check( nf90_put_att( this % ioVars % ncid, this % ioVars % xVarID, "long_name", "array of t-grid longitudes" ) )
      CALL Check( nf90_put_att( this % ioVars % ncid, this % ioVars % xVarID, "units", "degrees_east" ) )
      CALL Check( nf90_put_att( this % ioVars % ncid, this % ioVars % xVarID, "_FillValue", fillValue) )
      CALL Check( nf90_put_att( this % ioVars % ncid, this % ioVars % xVarID, "missing_value", fillValue) )


      ! Set up the tracer field names based on the model type
      IF( modelType == ImpulseField )THEN
      
         DO i = 1, this % nTracers
            WRITE( tracerid, '(I2.2)') i
            CALL Check( nf90_def_var( this % ioVars % ncid, "IRF_"//tracerid, NF90_BYTE,&
                                      (/ this % ioVars % xDimID, this % ioVars % yDimID, this % ioVars % zDimID /), &
                                      this % ioVars % tracerVarID(i) ) )
         ENDDO
         
      
      ELSE
 
         DO i = 1, this % nTracers
            WRITE( tracerid, '(I2.2)') this % tracerIds(i)-1
            CALL Check( nf90_def_var( this % ioVars % ncid, "DyeTracer_"//tracerid, NF90_FLOAT, &
                                      (/ this % ioVars % xDimID, this % ioVars % yDimID, this % ioVars % zDimID, this % ioVars % recDimID /), &
                                      this % ioVars % tracerVarID(i) ) )
   
            CALL Check( nf90_put_att( this % ioVars % ncid, this % ioVars % tracerVarID(i), "long_name", &
                                      "Dye tracer concentration of tracer "//tracerid ) )
            CALL Check( nf90_put_att( this % ioVars % ncid, this % ioVars % tracerVarID(i), "units", "%" ) )
            CALL Check( nf90_put_att( this % ioVars % ncid, this % ioVars % tracerVarID(i), "coordinates", &
                                      "TLONG TLAT z_t" ) )
            CALL Check( nf90_put_att( this % ioVars % ncid, this % ioVars % tracerVarID(i), "_FillValue", &
                                      fillValue) )
            CALL Check( nf90_put_att( this % ioVars % ncid, this % ioVars % tracerVarID(i), "missing_value", &
                                      fillValue) )

            CALL Check( nf90_def_var( this % ioVars % ncid, "Volume_"//tracerid, NF90_FLOAT, &
                                      (/ this % ioVars % xDimID, this % ioVars % yDimID, this % ioVars % zDimID, this % ioVars % recDimID /), &
                                       this % ioVars % volumeVarID(i) )  )

            CALL Check( nf90_put_att( this % ioVars % ncid, this % ioVars % volumeVarID(i), "long_name", &
                                      "Fractional change of fluid volume" ) )
            CALL Check( nf90_put_att( this % ioVars % ncid, this % ioVars % volumeVarID(i), "units", "" ) )
            CALL Check( nf90_put_att( this % ioVars % ncid, this % ioVars % volumeVarID(i), "coordinates", &
                                      "TLONG TLAT z_t" ) )
            CALL Check( nf90_put_att( this % ioVars % ncid, this % ioVars % volumeVarID(i), "_FillValue", &
                                      fillValue) )
            CALL Check( nf90_put_att( this % ioVars % ncid, this % ioVars % volumeVarID(i), "missing_value", &
                                      fillValue) )

         ENDDO

         IF( initOn ) THEN     

            DO i = 1, this % nTracers
               WRITE( tracerid, '(I2.2)') this % tracerIds(i)-1
               CALL Check( nf90_def_var( this % ioVars % ncid, "Source_"//tracerid, NF90_FLOAT, &
                                         (/ this % ioVars % xDimID, this % ioVars % yDimID, this % ioVars % zDimID /), &
                                         this % ioVars % sourceVarID(i) ) )
               CALL Check( nf90_def_var( this % ioVars % ncid, "rFac_"//tracerid, NF90_FLOAT, &
                                         (/ this % ioVars % xDimID, this % ioVars % yDimID, this % ioVars % zDimID /), &
                                         this % ioVars % rFacVarid(i) ) )
               CALL Check( nf90_def_var( this % ioVars % ncid, "mask_"//tracerid, NF90_FLOAT, &
                                         (/ this % ioVars % xDimID, this % ioVars % yDimID, this % ioVars % zDimID /), &
                                         this % ioVars % maskVarID(i) ) )
   
               CALL Check( nf90_put_att( this % ioVars % ncid, this % ioVars % sourceVarID(i), "long_name", &
                                         "Relaxation field of tracer "//tracerid ) )
               CALL Check( nf90_put_att( this % ioVars % ncid, this % ioVars % sourceVarID(i), "units", "" ) )
               CALL Check( nf90_put_att( this % ioVars % ncid, this % ioVars % sourceVarID(i), "coordinates", &
                                         "TLONG TLAT z_t" ) )
               CALL Check( nf90_put_att( this % ioVars % ncid, this % ioVars % sourceVarID(i), "_FillValue", &
                                         fillValue) )
               CALL Check( nf90_put_att( this % ioVars % ncid, this % ioVars % sourceVarID(i), "missing_value", &
                                         fillValue) )
   
               CALL Check( nf90_put_att( this % ioVars % ncid, this % ioVars % rFacVarid(i), "long_name", &
                                         "Relaxation frequency of tracer "//tracerid ) )
               CALL Check( nf90_put_att( this % ioVars % ncid, this % ioVars % rFacVarid(i), "units", "" ) )
               CALL Check( nf90_put_att( this % ioVars % ncid, this % ioVars % rFacVarid(i), "coordinates", &
                                         "TLONG TLAT z_t" ) )
               CALL Check( nf90_put_att( this % ioVars % ncid, this % ioVars % rFacVarid(i), "_FillValue", &
                                         fillValue) )
               CALL Check( nf90_put_att( this % ioVars % ncid, this % ioVars % rFacVarid(i), "missing_value", &
                                         fillValue) )
   
               CALL Check( nf90_put_att( this % ioVars % ncid, this % ioVars % maskVarID(i), "long_name", &
                                         "Source-region mask of tracer "//tracerid ) )
               CALL Check( nf90_put_att( this % ioVars % ncid, this % ioVars % maskVarID(i), "units", "" ) )
               CALL Check( nf90_put_att( this % ioVars % ncid, this % ioVars % maskVarID(i), "coordinates", &
                                         "TLONG TLAT z_t" ) )
               CALL Check( nf90_put_att( this % ioVars % ncid, this % ioVars % maskVarID(i), "_FillValue", &
                                         fillValue) )
               CALL Check( nf90_put_att( this % ioVars % ncid, this % ioVars % maskVarID(i), "missing_value", &
                                         fillValue) )
            ENDDO

         ENDIF              
      
      ENDIF

      ! End the Define Mode
      CALL Check( nf90_enddef(this % ioVars % ncid) )
      ! 
      CALL Check( nf90_put_var( this % ioVars % ncid, this % ioVars % zVarID, mesh % z ) )
      CALL Check( nf90_put_var( this % ioVars % ncid, this % ioVars % yVarID, mesh % tLat ) )
      CALL Check( nf90_put_var( this % ioVars % ncid, this % ioVars % xVarID, mesh % tLon ) )

      INFO('Done.')

 END SUBROUTINE InitializeForNetCDFWrite_POP_Native
! 
 SUBROUTINE LoadNetCDF_POP_IRF( this, filename )
   IMPLICIT NONE
   CLASS( POP_IRF ), INTENT(inout) :: this
   CHARACTER(*), INTENT(in)        :: filename
   ! Local
   INTEGER :: i, ncid
   INTEGER :: start(1:3), recCount(1:3)
   CHARACTER(2) :: tracerid
   INTEGER :: varid(1:this % nIRF+1)

      start    = (/1, 1, 1/)
      recCount = (/this % nX, this % nY, this % nZ/)

      CALL Check( nf90_open( TRIM(filename), nf90_nowrite, ncid) )

      
      DO i = 1, this % nIRF
         WRITE( tracerid, '(I2.2)') i
         CALL Check( nf90_inq_varid( ncid, "ADV_3D_IRF_"//tracerid, &
                                     varID(i) ) )
      ENDDO
      CALL Check( nf90_inq_varid( ncid, "VDC_S", &
                                  varID(this % nIRF+1) ) )


      DO i = 1, this % nIRF+1
         CALL Check( nf90_get_var( ncid, &
                                   varID(i), &
                                   this % irf(:,:,:,i), &
                                   start, recCount ) )
      ENDDO

      CALL Check( nf90_close( ncid ) )

 END SUBROUTINE LoadNetCDF_POP_IRF

 SUBROUTINE InitializeForNetCDFRead_POP_Native( this, modelType, filename, initOn )
#undef __FUNC__
#define __FUNC__ "InitializeForNetCDFRead_POP_Native"
   IMPLICIT NONE
   CLASS( POP_Native ), INTENT(inout) :: this
   INTEGER, INTENT(in)             :: modelType
   CHARACTER(*), INTENT(in)        :: filename
   LOGICAL, INTENT(in)             :: initOn
   ! Local
   INTEGER :: i
   CHARACTER(2) :: tracerid

      INFO('Initializing NetCDF file for read. '//TRIM(filename))
      ! Create the netcdf file and generate a file handle referenced by the
      ! integer "this % ioVars % ncid"
      CALL Check( nf90_open( TRIM(filename), nf90_nowrite, this % ioVars % ncid ) )
      ! Create variables -- here we need to create arrays for the dimensions
      CALL Check( nf90_inq_varid( this % ioVars % ncid, "z_t", this % ioVars % zVarID ) )
      CALL Check( nf90_inq_varid( this % ioVars % ncid, "TLAT", this % ioVars % yVarID ) )
      CALL Check( nf90_inq_varid( this % ioVars % ncid, "TLONG", this % ioVars % xVarID ) )

      ! Set up the tracer field names based on the model type
      IF( modelType == ImpulseResponseField )THEN
      
         DO i = 1, this % nTracers-1
            WRITE( tracerid, '(I2.2)') i
            CALL Check( nf90_inq_varid( this % ioVars % ncid, "ADV_3D_IRF_"//tracerid, &
                                        this % ioVars % tracerVarID(i) ) )
         ENDDO
         CALL Check( nf90_inq_varid( this % ioVars % ncid, "VDC_S", &
                                     this % ioVars % tracerVarID(this % nTracers) ) )

      ELSEIF( modelType == ImpulseField )THEN
      
         DO i = 1, this % nTracers
            WRITE( tracerid, '(I2.2)')i
            CALL Check( nf90_inq_varid( this % ioVars % ncid, "ADV_3D_IRF_"//tracerid, &
                                        this % ioVars % tracerVarID(i) ) )
         ENDDO
                             
      
      ELSE

         IF( initOn )THEN
      
            DO i = 1, this % nTracers
               WRITE( tracerid, '(I2.2)') this % tracerIds(i)-1
               INFO('Obtain variable ID for DyeTracer_'//TRIM(tracerid))
               CALL Check( nf90_inq_varid( this % ioVars % ncid, "DyeTracer_"//tracerid, &
                                         this % ioVars % tracerVarID(i) ) )
               INFO('Obtain variable ID for Volume_'//TRIM(tracerid))
               CALL Check( nf90_inq_varid( this % ioVars % ncid,"Volume_"//tracerid, &
                                           this % ioVars % volumeVarID(i) ) )
               INFO('Obtain variable ID for Source_'//TRIM(tracerid))
               CALL Check( nf90_inq_varid( this % ioVars % ncid, "Source_"//tracerid, &
                                         this % ioVars % sourceVarID(i) ) )
               INFO('Obtain variable ID for rFac_'//TRIM(tracerid))
               CALL Check( nf90_inq_varid( this % ioVars % ncid, "rFac_"//tracerid, &
                                         this % ioVars % rFacVarid(i) ) )
               INFO('Obtain variable ID for mask_'//TRIM(tracerid))
               CALL Check( nf90_inq_varid( this % ioVars % ncid, "mask_"//tracerid, &
                                         this % ioVars % maskVarID(i) ) )
            ENDDO

         ELSE

            DO i = 1, this % nTracers
               WRITE( tracerid, '(I2.2)') this % tracerIds(i)-1
               INFO('Obtain variable ID for DyeTracer_'//TRIM(tracerid))
               CALL Check( nf90_inq_varid( this % ioVars % ncid, "DyeTracer_"//tracerid, &
                                         this % ioVars % tracerVarID(i) ) )
               INFO('Obtain variable ID for Volume_'//TRIM(tracerid))
               CALL Check( nf90_inq_varid( this % ioVars % ncid,"Volume_"//tracerid, &
                                           this % ioVars % volumeVarID(i) ) )
            ENDDO
         
         ENDIF 

      
      ENDIF

      INFO('Done.')

 END SUBROUTINE InitializeForNetCDFRead_POP_Native
!   
       
 SUBROUTINE FinalizeNetCDF_POP_Native( this )
    IMPLICIT NONE
    CLASS( POP_Native ) :: this
 
       CALL Check( nf90_close( this % ioVars % ncid ) )
       
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

      CALL Check( nf90_inq_varid( ncid, "SSH",varid ) )
      CALL Check( nf90_get_var( ncid, &
                                varid, &
                                this % volume(:,:,1,1), &
                                start2D, recCount2D ) )

!      CALL Check( nf90_inq_varid( ncid, "TEMP",varid ) )
!      CALL Check( nf90_get_var( ncid, &
!                                varid, &
!                                this % temperature, &
!                                start3D, recCount3D ) )
!
!      CALL Check( nf90_inq_varid( ncid, "SALT",varid ) )
!      CALL Check( nf90_get_var( ncid, &
!                                varid, &
!                                this % salinity, &
!                                start3D, recCount3D ) )
!
!      CALL Check( nf90_inq_varid( ncid, "PD",varid ) )
!      CALL Check( nf90_get_var( ncid, &
!                                varid, &
!                                this % density, &
!                                start3D, recCount3D ) )
!

      CALL Check( nf90_close( ncid ) )

 END SUBROUTINE LoadOceanState_POP_Native
!
 SUBROUTINE WriteOceanState_POP_Native( this, mesh, filename )
   IMPLICIT NONE
   CLASS( POP_Native ), INTENT(inout) :: this
   TYPE( POP_Mesh ), INTENT(in)       :: mesh
   CHARACTER(*), INTENT(in)           :: filename
   ! Local
   INTEGER :: ncid, zDimID, xDimID, yDimID
   INTEGER :: varid_ssh, varid_temp, varid_salt, varid_pd
   INTEGER :: start2D(1:2), recCount2D(1:2)
   INTEGER :: start3D(1:3), recCount3D(1:3)

      CALL Check( nf90_create( PATH=TRIM(filename),&
                               CMODE=OR(nf90_clobber,nf90_64bit_offset),&
                               NCID=ncid ) )
      ! Create the dimensions - the dimension names are currently chosen based
      CALL Check( nf90_def_dim( ncid, "z_t", mesh % nZ, zDimID ) ) 
      CALL Check( nf90_def_dim( ncid, "nlon", mesh % nX, xDimID ) ) 
      CALL Check( nf90_def_dim( ncid, "nlat", mesh % nY, yDimID ) ) 


      start2D    = (/1, 1/)
      recCount2D = (/mesh % nX, mesh % nY/)
      start3D    = (/1, 1, 1/)
      recCount3D = (/mesh % nX, mesh % nY, mesh % nZ/)

!      CALL Check( nf90_def_var( ncid, "SSH", NF90_FLOAT,&
!                                (/ xDimID, yDimID /), &
!                                 varid_ssh ) )
!
!      CALL Check( nf90_def_var( ncid, "TEMP", NF90_FLOAT,&
!                                (/ xDimID, yDimID, zDimID /), &
!                                 varid_temp ) )
!
!      CALL Check( nf90_def_var( ncid, "SALT", NF90_FLOAT,&
!                                (/ xDimID, yDimID, zDimID /), &
!                                 varid_salt ) )
!
!      CALL Check( nf90_def_var( ncid, "PD", NF90_FLOAT,&
!                                (/xDimID, yDimID, zDimID/),&
!                                 varid_pd ) )
!

!      CALL Check( nf90_enddef(ncid) )
!
!      CALL Check( nf90_put_var( ncid, &
!                                varid_ssh, &
!                                this % volume(:,:,1), &
!                                start2D, recCount2D ) )
!
!      CALL Check( nf90_put_var( ncid, &
!                                varid_temp, &
!                                this % temperature, &
!                                start3D, recCount3D ) )
!
!      CALL Check( nf90_put_var( ncid, &
!                                varid_salt, &
!                                this % salinity, &
!                                start3D, recCount3D ) )
!
!      CALL Check( nf90_put_var( ncid, &
!                                varid_pd, &
!                                this % density, &
!                                start3D, recCount3D ) )
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
         CALL Check( nf90_put_var( this % ioVars % ncid, &
                                   this % ioVars % sourceVarID(i), &
                                   this % source(:,:,:,i), &
                                   start, recCount ) )      
         CALL Check( nf90_put_var( this % ioVars % ncid, &
                                   this % ioVars % rfacVarID(i), &
                                   this % rfac(:,:,:,i), &
                                   start, recCount ) )      
         CALL Check( nf90_put_var( this % ioVars % ncid, &
                                   this % ioVars % maskVarID(i), &
                                   this % mask(:,:,:,i), &
                                   start, recCount ) )      
      ENDDO

 END SUBROUTINE WriteSourceEtcNETCDF_POP_Native
!
 SUBROUTINE WriteNetCDFRecord_POP_Native( this, mesh, recordID )
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
         CALL Check( nf90_put_var( this % ioVars % ncid, &
                                   this % ioVars % tracerVarID(i), &
                                   this % tracer(:,:,:,i), &
                                   start, recCount ) )      
         CALL Check( nf90_put_var( this % ioVars % ncid, &
                                   this % ioVars % volumeVarID(i), &
                                   this % volume(:,:,:,i), &
                                   start, recCount ) )      
      ENDDO

 END SUBROUTINE WriteNetCDFRecord_POP_Native
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
            CALL Check( nf90_get_var( this % ioVars % ncid, &
                                      this % ioVars % sourceVarID(i), &
                                      this % source(:,:,:,i), &
                                      start, recCount ) )
            CALL Check( nf90_get_var( this % ioVars % ncid, &
                                      this % ioVars % rfacVarID(i), &
                                      this % rfac(:,:,:,i), &
                                      start, recCount ) )
            CALL Check( nf90_get_var( this % ioVars % ncid, &
                                      this % ioVars % maskVarID(i), &
                                      this % mask(:,:,:,i), &
                                      start, recCount ) )
         ENDDO

 END SUBROUTINE ReadSourceEtcNETCDF_POP_Native
!
 SUBROUTINE ReadNetCDFRecord_POP_Native( this, mesh, recordID )
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
            CALL Check( nf90_get_var( this % ioVars % ncid, &
                                      this % ioVars % tracerVarID(i), &
                                      this % tracer(:,:,:,i), &
                                      start, recCount ) )
            CALL Check( nf90_get_var( this % ioVars % ncid, &
                                      this % ioVars % volumeVarID(i), &
                                      this % volume(:,:,:,i), &
                                      start, recCount ) )
         ENDDO

 END SUBROUTINE ReadNetCDFRecord_POP_Native
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
            CALL Check( nf90_get_var( this % ioVars % ncid, &
                                      this % ioVars % tracerVarID(i), &
                                      this % tracer(:,:,:,i), &
                                      start, recCount ) )
         ENDDO

 END SUBROUTINE LoadTracerFromNetCDF_POP_Native
!
 SUBROUTINE WriteTracerToNetCDF_POP_Native( this, mesh )
   IMPLICIT NONE
   CLASS( POP_Native ), INTENT(inout) :: this
   TYPE( POP_Mesh ), INTENT(in)       :: mesh
   ! Local
   INTEGER :: start(1:3), recCount(1:3)
   INTEGER :: i

         start    = (/1, 1, 1/)
         recCount = (/mesh % nX, mesh % nY, mesh % nZ/)

         DO i = 1, this % nTracers
            CALL Check( nf90_put_var( this % ioVars % ncid, &
                                      this % ioVars % tracerVarID(i), &
                                      this % tracer(:,:,:,i), &
                                      start, recCount ) )
         ENDDO

 END SUBROUTINE WriteTracerToNetCDF_POP_Native

END MODULE POP_Native_Class
