PROGRAM FEOTSInitialize


! src/POP/
USE POP_FEOTS_Class
USE POP_Native_Class
USE POP_Params_Class


IMPLICIT NONE

   TYPE( POP_FEOTS ) :: feots


      CALL feots % Build( )

      CALL InitialConditions( feots )

      !  //////////////////////////////////////////// File I/O  //////////////////////////////////////////////////////// !
      CALL feots % nativeSol % InitializeForNetCDFWrite( feots % params % TracerModel, &
                                                         feots % mesh, &
                                                         TRIM(feots % params % outputDirectory)//'Tracer.init.nc' )
      CALL feots % nativeSol % WriteNetCDFRecord( feots % mesh, 1 )
      CALL feots % nativeSol % WriteSourceEtcNetCDF( feots % mesh )

      CALL feots % nativeSol % FinalizeNetCDF( )
      ! //////////////////////////////////////////////////////////////////////////////////////////////////////////////// !

      CALL feots % Trash( )

CONTAINS

 SUBROUTINE InitialConditions( myfeots )
 ! Sets the initial tracer distributions
   IMPLICIT NONE
   TYPE( POP_FEOTS ), INTENT(inout) :: myFeots
   ! Local
   INTEGER  :: i, j, k, iTracer, m, dof
   REAL(prec) :: x, y, z, rfMax, s, T, rho


      rfMax = 1.0_prec/43200.0_prec
      DO k = 1, myFeots % mesh % nZ  

         z = myFeots % mesh % z(k)

         DO j = 1, myFeots % mesh % nY
            DO i = 1, myFeots % mesh % nX 
 
               x = myFeots % mesh % tLon(i,j)
               y = myFeots % mesh % tLat(i,j)

               ! Agulhas
               myFeots % nativeSol % tracer(i,j,k,:)   = 0.0_prec 

               myFeots % nativeSol % source(i,j,k,:)   = 0.0_prec !5.0_prec
               myFeots % nativeSol % rFac(i,j,k,:)     = 0.0_prec !rFMax*(exp( -0.5_prec*( (x-31.5_prec)**2 +(y+31.0_prec)**2 )/(0.25_prec) ))
      
            ENDDO
         ENDDO
      ENDDO

      ! Set any prescribed cells here
      DO m = 1, myFeots % regionalMaps % nPCells

         dof = myFeots % regionalMaps % prescribedCells(m)
         i   = myFeots % regionalMaps % dofToLocalIJK(1,dof)
         j   = myFeots % regionalMaps % dofToLocalIJK(2,dof)
         k   = myFeots % regionalMaps % dofToLocalIJK(3,dof)

         x = myFeots % mesh % tLon(i,j)
         y = myFeots % mesh % tLat(i,j)
         myFEOTS % nativeSol % tracer(i,j,k,:) = 0.0_prec

      ENDDO

 END SUBROUTINE InitialConditions
!
END PROGRAM FEOTSInitialize
