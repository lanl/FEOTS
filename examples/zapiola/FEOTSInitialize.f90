PROGRAM FEOTSInitialize


! src/POP/
USE POP_FEOTS_Class
USE POP_Native_Class
USE POP_Params_Class


IMPLICIT NONE

   TYPE( POP_FEOTS ) :: feots
   INTEGER :: mpiErr, myRank, nProcs

#ifdef HAVE_MPI
      CALL MPI_INIT( mpiErr )
      CALL MPI_COMM_RANK( MPI_COMM_WORLD, myRank, mpiErr )
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD, nProcs, mpiErr )
#else
      myRank = 0
      nProcs = 1
#endif

      CALL feots % Build( myRank, nProcs )

      print*, 'Check!'
      
      CALL SourceTerms( feots )

      !  //////////////////////////////////////////// File I/O  //////////////////////////////////////////////////////// !
      IF( myRank == 0 )THEN
        CALL feots % nativeSol % InitializeForNetCDFWrite( feots % params % TracerModel, &
                                                           feots % mesh, &
                                                           TRIM(feots % params % outputDirectory)//'Tracer.init.nc', &
                                                           .TRUE. )
        CALL feots % nativeSol % WriteNetCDFRecord( feots % mesh, 1 )
        CALL feots % nativeSol % WriteSourceEtcNetCDF( feots % mesh )
  
        CALL feots % nativeSol % FinalizeNetCDF( )
      ENDIF
      ! //////////////////////////////////////////////////////////////////////////////////////////////////////////////// !

      CALL feots % Trash( )

#ifdef HAVE_MPI
      CALL MPI_FINALIZE( mpiErr )
#endif

CONTAINS

 SUBROUTINE SourceTerms( myFeots )
 ! Sets the source terms and the "relaxation factor" for each tracer
 !
   IMPLICIT NONE
   TYPE( POP_FEOTS ), INTENT(inout) :: myFeots
   ! Local
   INTEGER  :: i, j, k, m, nn, iTracer

      IF ( myRank == 0 )THEN

      DO iTracer = 1, myFeots % params % nTracers
        PRINT*, 'RANK, nPCELLS', iTracer, myFeots % regionalMaps % bMap(iTracer) % nPCells
        DO m = 1, myFeots % regionalMaps % bMap(iTracer) % nPCells

               nn = myFeots % regionalMaps % bMap(iTracer) % prescribedCells(m)
               i   = myFeots % regionalMaps % dofToLocalIJK(1,nn)
               j   = myFeots % regionalMaps % dofToLocalIJK(2,nn)
               k   = myFeots % regionalMaps % dofToLocalIJK(3,nn)

! Southern boundary

         IF( k .le. 63 .AND. iTracer==1 ) THEN

             print*, 'RANK, I,J,K', iTracer, i, j, k
            myFeots % nativeSol % tracer(i,j,k,iTracer) = 1.0_prec

         ELSEIF ( k .gt. 63 .AND. iTracer==2 ) THEN

             print*, 'RANK, I,J,K', iTracer, i, j, k
            myFeots % nativeSol % tracer(i,j,k,iTracer) = 1.0_prec

         ENDIF

! Eastern boundary

         IF ( k .le. 63 .AND. iTracer == 3 ) THEN

             print*, 'RANK, I,J,K', iTracer, i, j, k
            myFeots % nativeSol % tracer(i,j,k,iTracer) = 1.0_prec

         ELSEIF ( k .gt. 63 .AND.  iTracer == 4 ) THEN

             print*, 'RANK, I,J,K', iTracer, i, j, k
            myFeots % nativeSol % tracer(i,j,k,iTracer) = 1.0_prec

         ENDIF

! Northern boundary

         IF (k .le. 63 .AND. iTracer == 5 ) THEN

             print*, 'RANK, I,J,K', iTracer, i, j, k
            myFeots % nativeSol % tracer(i,j,k,iTracer) = 1.0_prec

         ELSEIF ( k .gt. 63 .AND. iTracer == 6 ) THEN

             print*, 'RANK, I,J,K', iTracer, i, j, k
            myFeots % nativeSol % tracer(i,j,k,iTracer) = 1.0_prec

         ENDIF


        ENDDO
    ENDDO
    ENDIF
 END SUBROUTINE SourceTerms
!  

END PROGRAM FEOTSInitialize
