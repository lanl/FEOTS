PROGRAM FEOTSInitialize


! src/POP/
USE POP_FEOTS_Class
USE POP_Native_Class
USE POP_Params_Class


IMPLICIT NONE

   TYPE( POP_FEOTS ) :: feots


      CALL feots % Build( myRank = 0, nProcs = 1 )

      CALL InitialConditions( feots )

      !  //////////////////////////////////////////// File I/O  //////////////////////////////////////////////////////// !
      CALL feots % nativeSol % InitializeForNetCDFWrite( feots % params % TracerModel, &
                                                         feots % mesh, &
                                                         TRIM(feots % params % outputDirectory)//'Tracer.init.nc', &
                                                         .TRUE. )
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
   INTEGER  :: i, j, k, iTracer, m, dof, iMask, iLayer
   REAL(prec) :: x, y, z, rfMax,trackingVar


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
               myFeots % nativeSol % rFac(i,j,k,:)     = 0.0_prec !rFMax*(exp(
! -0.5_prec*( (x-31.5_prec)**2 +(y+31.0_prec)**2 )/(0.25_prec) ))

            ENDDO
         ENDDO
      ENDDO

      ! Set any prescribed cells here
      DO iMask = 1, myFeots % regionalMaps % nMasks
         DO iLayer = 1, myFeots % params % nLayers
            iTracer = iLayer + (iMask-1)*(myFeots % params % nLayers)

            DO m = 1, myFeots % regionalMaps % bMap(iMask) % nPCells

               dof = myFeots % regionalMaps % bMap(iMask) % prescribedCells(m)
               i   = myFeots % regionalMaps % dofToLocalIJK(1,dof)
               j   = myFeots % regionalMaps % dofToLocalIJK(2,dof)
               k   = myFeots % regionalMaps % dofToLocalIJK(3,dof)

               ! x = myFeots % mesh % tLon(i,j)
               ! y = myFeots % mesh % tLat(i,j)
               ! myFEOTS % nativeSol % tracer(i,j,k,iTracer) = 1.0_prec

               ! --------------------------------------------------------------------
               !***** This section of code provides an example of the water
               !***** mass tagging implementation. Remove this section if you're
               !***** not doing water mass tagging. Note that when water mass
               !***** tagging is enabled, this hard-setting is done during
               !forward
               !***** integration; assigning the hard set values in the initial
               !***** condition is redundant but helpful for consistency
               !checking.
               !
               ! Assign the tracer value according to the water mass layer and
               ! the
               ! appropriate boundary mask
               trackingVar = myFeots % nativeSol % temperature(i,j,k)*myFeots % statemask(1) +&
                             myFeots % nativeSol % salinity(i,j,k)*myFeots % statemask(2) +&
                             myFeots % nativeSol % density(i,j,k)*myFeots % statemask(3)

               IF( trackingVar >= myFeots % stateLowerBound(iLayer) .AND. &
                   trackingVar < myFeots % stateUpperBound(iLayer) )THEN
                   myFeots % nativeSol % tracer(i,j,k,iTracer) = 1.0_prec
               ENDIF
               ! -----------------------------------------------------------------
            ENDDO
          ENDDO
      ENDDO

 END SUBROUTINE InitialConditions

END PROGRAM FEOTSInitialize
