! GenMask.f90
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

PROGRAM GenMask

! src/common/
USE CommonRoutines
! src/POP/
USE POP_Params_Class
USE POP_Mesh_Class
USE FEOTS_CLI_Class


IMPLICIT NONE
   REAL(prec), PARAMETER :: prescribed_width = 0.10_prec

   TYPE( POP_Params )   :: params
   TYPE( POP_Mesh )     :: mesh
   INTEGER              :: i, j, iMask, fUnit
   INTEGER, ALLOCATABLE :: maskField(:,:,:)
   INTEGER, ALLOCATABLE :: regionMask(:,:)
   CHARACTER(400)       :: ncfile
   CHARACTER(200)       :: IRFfile
   REAL(prec)           :: x, y, r
   INTEGER              :: nMasks
   TYPE(FEOTS_CLI) :: cliParams

      CALL cliParams % GetCLIConf( )

      CALL params % Build( cliParams % paramFile )

      CALL mesh % Load( TRIM(cliParams % dbRoot)//'/mesh/mesh.nc'  )
  
      nMasks = params % nTracers

      ALLOCATE( maskField(1:mesh % nX,1:mesh % nY,1:nMasks) )
      ALLOCATE( regionMask(1:mesh % nX,1:mesh % nY) )

      maskField  = 0
      regionMask = 0

!      CALL LoadRegionMask( mesh, regionMask, TRIM(cliParams % irfFile) )


         DO j = 1, mesh % nY
            DO i = 1, mesh % nX

               x = mesh % tLon(i,j)
               y = mesh % tLat(i,j)
             
               IF( x >= 180.0_prec )THEN
                  x = x -360.0_prec
               ENDIF
              
               IF( y >= params % south .AND. y <= params % north .AND. &
                   x >= params % west .AND. x <= params % east)THEN
 
                 ! Set interior points for the mesh
                 maskfield(i,j,1:6) = 1

                 IF( y < params % south + prescribed_width )THEN
                   maskfield(i,j,1:2) = -1 ! prescribed
                 ENDIF

                 IF( x > params % east - prescribed_width )THEN
                   maskfield(i,j,3:4) = -1 ! prescribed
                 ENDIF

                 IF( y > params % north - prescribed_width )THEN
                   maskfield(i,j,5:6) = -1 ! prescribed
                 ENDIF
              
               ENDIF

            ENDDO
         ENDDO

      
      CALL WriteMaskField( mesh, maskField, TRIM(cliParams % outdir)//'/mask.nc' )
      CALL mesh % Trash( )
      DEALLOCATE( regionMask, maskField )

CONTAINS

 SUBROUTINE WriteMaskField( mesh, maskfield, maskfile )

   IMPLICIT NONE
   TYPE( POP_Mesh ), INTENT(inout)      :: mesh
   INTEGER, INTENT(in)                  :: maskfield(1:mesh % nX, 1:mesh %nY,1:nMasks)
   CHARACTER(*), INTENT(in)             :: maskfile
   ! Local
   INTEGER :: start(1:2), recCount(1:2)
   INTEGER :: ncid, varid(1:nMasks), nMaskid, x_dimid, y_dimid, iMask
   CHARACTER(3) :: maskChar

      PRINT*, 'Writing mask to '//TRIM(maskfile)
      start    = (/1, 1/)
      recCount = (/mesh % nX, mesh % nY/)

      CALL Check( nf90_create( PATH=TRIM(maskfile),&
                               CMODE=OR(nf90_clobber,nf90_64bit_offset),&
                               NCID=ncid ) )
      CALL Check( nf90_def_dim( ncid, "nlon", mesh % nX, x_dimid ) )
      CALL Check( nf90_def_dim( ncid, "nlat", mesh % nY, y_dimid ) )
      CALL Check( nf90_def_var( ncid, "nMasks",NF90_INT, nMaskid ) )

      DO iMask = 1, nMasks
         WRITE( maskChar,'(I3.3)' )iMask
         CALL Check( nf90_def_var( ncid, "mask"//maskChar,NF90_INT,&
                               (/ x_dimid, y_dimid /),&
                                varid(iMask) ) )

         CALL Check( nf90_put_att( ncid, varid(iMask), "long_name", "Domain Mask" ) )
         CALL Check( nf90_put_att( ncid, varid(iMask), "units", "" ) )
      ENDDO

      CALL Check( nf90_enddef(ncid) )
 
      CALL Check( nf90_put_var(ncid, nMaskid, nMasks ) ) 

      DO iMask = 1, nMasks
         CALL Check( nf90_put_var( ncid, &
                                varid(iMask), &
                                maskfield(:,:,iMask), &
                                start, recCount ) )
      ENDDO

      CALL Check( nf90_close( ncid ) )


 END SUBROUTINE WriteMaskField
SUBROUTINE LoadRegionMask(  mesh, regionmask, filename )
   IMPLICIT NONE
   TYPE( POP_Mesh ), INTENT(in)       :: mesh
   INTEGER, INTENT(out)               :: regionmask(1:mesh % nX, 1:mesh % nY)
   CHARACTER(*), INTENT(in)           :: filename
   ! Local
   INTEGER :: ncid, varid
   INTEGER :: start2D(1:2), recCount2D(1:2)

      CALL Check( nf90_open( TRIM(filename), nf90_nowrite, ncid ) )

      start2D    = (/1, 1/)
      recCount2D = (/mesh % nX, mesh % nY/)

      PRINT*, 'Loading Region Mask'
      CALL Check( nf90_inq_varid( ncid, "REGION_MASK",varid ) )
      CALL Check( nf90_get_var( ncid, &
                                varid, &
                                regionMask, &
                                start2D, recCount2D ) )

      PRINT*, 'DONE'

      CALL Check( nf90_close( ncid ) )

 END SUBROUTINE LoadRegionMask


END PROGRAM GenMask

