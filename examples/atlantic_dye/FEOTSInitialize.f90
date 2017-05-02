PROGRAM FEOTSInitialize


! src/POP/
USE POP_FEOTS_Class
USE POP_Native_Class
USE POP_Params_Class


IMPLICIT NONE

   TYPE( POP_FEOTS ) :: feots


      CALL feots % Build( )

      CALL InitialConditions( feots )

      CALL SourceTerms( feots )

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
   INTEGER  :: i, j, k
   REAL(prec) :: x, y, z


      DO k = 1, myFeots % mesh % nZ  

         z = myFeots % mesh % z(k)

         DO j = 1, myFeots % mesh % nY
            DO i = 1, myFeots % mesh % nX 
 
               x = myFeots % mesh % tLon(i,j)
               y = myFeots % mesh % tLat(i,j)

               ! Agulhas
               myFeots % nativeSol % tracer(i,j,k,:) = 0.0_prec ! 5.0_prec*exp( -0.5_prec*((x-31.0_prec)**2 +(y+32.0_prec)**2 )/(0.25_prec) )

            ENDDO
         ENDDO
      ENDDO

 END SUBROUTINE InitialConditions
!
 SUBROUTINE SourceTerms( myFeots )
 ! Sets the source terms and the "relaxation factor" for each tracer
 !
   IMPLICIT NONE
   TYPE( POP_FEOTS ), INTENT(inout) :: myFeots
   ! Local
   INTEGER  :: i, j, k, m, nn
   REAL(prec) :: rFmax, x, y, z

      rFmax = 1.0_prec/86400.0_prec

      DO k = 1, myFeots % mesh % nZ  

         z = myFeots % mesh % z(k)

         DO j = 1, myFeots % mesh % nY
            DO i = 1, myFeots % mesh % nX 
 
               x = myFeots % mesh % tLon(i,j)
               y = myFeots % mesh % tLat(i,j)

              !  ! Agulhas
              !  myFeots % nativeSol % source(i,j,k,:) = 5.0_prec*exp( -0.5_prec*( (x-31.5_prec)**2 +(y+31.0_prec)**2 )/(0.25_prec) )
              !  myFeots % nativeSol % rFac(i,j,k,:)   = rFMax*(exp( -0.5_prec*( (x-31.5_prec)**2 +(y+31.0_prec)**2 )/(0.25_prec) ))**2

               ! Agulhas
               myFeots % nativeSol % source(i,j,k,:) = 0.0_prec
               myFeots % nativeSol % rFac(i,j,k,:)   = 0.0_prec
            ENDDO
         ENDDO
      ENDDO


      DO m = 1, myFeots % regionalMaps % nBCells

         nn = myFeots % regionalMaps % boundaryCells(m)

         i = myFeots % regionalMaps % dofToLocalIJK (1,nn)      ! lon
         j = myFeots % regionalMaps % dofToLocalIJK (2,nn)      ! lat
         k = myFeots % regionalMaps % dofToLocalIJK (3,nn)      ! depth

         x = myFeots % mesh % tLon(i,j)
         y = myFeots % mesh % tLat(i,j)


         IF (x .gt. 29.0_prec .and. y .lt. 0.0_prec) THEN
            
            myFeots % nativeSol % tracer(i,j,k,:) = 1.0_prec

         ENDIF

      ENDDO
 END SUBROUTINE SourceTerms
!  

END PROGRAM FEOTSInitialize
