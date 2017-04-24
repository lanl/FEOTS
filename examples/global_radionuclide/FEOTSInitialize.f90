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

      feots % nativeSol % tracer = (1.0_prec-feots % nativeSol % mask)*feots % nativeSol % hardset

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

               myFeots % nativeSol % tracer(i,j,k,:) = 0.0_prec 
               
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
   INTEGER  :: i, j, k
   REAL(prec) :: rFmax, x, y, z

      rFmax = 1.0_prec/86400.0_prec

      DO k = 1, myFeots % mesh % nZ  

         z = myFeots % mesh % z(k)

         DO j = 1, myFeots % mesh % nY
            DO i = 1, myFeots % mesh % nX 
 
               x = myFeots % mesh % tLon(i,j)
               y = myFeots % mesh % tLat(i,j)

               ! Particulate source
               myFeots % nativeSol % source(i,j,k,1) = 0.0_prec
               ! Radionuclide uniform source
               myFeots % nativeSol % source(i,j,k,2) = 1.0_prec*10.0_prec**(-6)
               myFeots % nativeSol % rFac(i,j,k,:)   = 0.0_prec 

               ! Set the Surface Tracer distribution
               IF( k == 1 )THEN
                  myFeots % nativeSol % mask(i,j,k,1)    = 0.0_prec
                  myFeots % nativeSol % hardSet(i,j,k,1) = 1.0_prec
               ENDIF


            ENDDO
         ENDDO
      ENDDO
 END SUBROUTINE SourceTerms
!  

END PROGRAM FEOTSInitialize
