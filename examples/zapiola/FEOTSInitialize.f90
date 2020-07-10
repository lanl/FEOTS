PROGRAM FEOTSInitialize


! src/POP/
USE POP_FEOTS_Class
USE POP_Native_Class
USE POP_Params_Class


IMPLICIT NONE

   TYPE( POP_FEOTS ) :: feots
   TYPE(FEOTS_CLI) :: cliParams
   INTEGER :: mpiErr, myRank, nProcs
   CHARACTER(5) :: rankChar

      CALL MPI_INIT( mpiErr )
      CALL MPI_COMM_RANK( MPI_COMM_WORLD, myRank, mpiErr )
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD, nProcs, mpiErr )

      WRITE( rankChar, '(I5.5)' ) myRank
      CALL cliParams % GetCLIConf( )
      CALL feots % Build( cliParams, myRank, nProcs )

      print*, 'Check!'
      
      CALL SourceTerms( feots )

      !  //////////////////////////////////////////// File I/O  //////////////////////////////////////////////////////// !
      CALL feots % nativeSol % InitializeForNetCDFWrite( feots % params % TracerModel, &
                                                         feots % mesh, &
                                                         TRIM(cliParams % outDir)//'/Tracer.'//rankChar//'.init.nc', &
                                                         .TRUE. )
      CALL feots % nativeSol % WriteNetCDFRecord( feots % mesh, 1 )
      CALL feots % nativeSol % WriteSourceEtcNetCDF( feots % mesh )
  
      CALL feots % nativeSol % FinalizeNetCDF( )
      ! //////////////////////////////////////////////////////////////////////////////////////////////////////////////// !

      CALL feots % Trash( )

      CALL MPI_FINALIZE( mpiErr )

CONTAINS

 SUBROUTINE SourceTerms( myFeots )
 ! Sets the source terms and the "relaxation factor" for each tracer
 !
   IMPLICIT NONE
   TYPE( POP_FEOTS ), INTENT(inout) :: myFeots
   ! Local
   INTEGER  :: i, j, k, m, nn, iTracer

      iTracer = 1  !1, myFeots % params % nTracers
        PRINT*, 'RANK, nPCELLS', iTracer, myFeots % regionalMaps % bMap(iTracer) % nPCells
        DO m = 1, myFeots % regionalMaps % bMap(iTracer) % nPCells

               nn = myFeots % regionalMaps % bMap(iTracer) % prescribedCells(m)
               i   = myFeots % regionalMaps % dofToLocalIJK(1,nn)
               j   = myFeots % regionalMaps % dofToLocalIJK(2,nn)
               k   = myFeots % regionalMaps % dofToLocalIJK(3,nn)

! Southern boundary

         IF( k .le. 63 .AND. myRank==1 ) THEN

             print*, 'RANK, I,J,K', iTracer, i, j, k
            myFeots % nativeSol % tracer(i,j,k,iTracer) = 1.0_prec

         ELSEIF ( k .gt. 63 .AND. myRank==2 ) THEN

             print*, 'RANK, I,J,K', iTracer, i, j, k
            myFeots % nativeSol % tracer(i,j,k,iTracer) = 1.0_prec

         ENDIF

! Eastern boundary

         IF ( k .le. 63 .AND. myRank == 3 ) THEN

             print*, 'RANK, I,J,K', iTracer, i, j, k
            myFeots % nativeSol % tracer(i,j,k,iTracer) = 1.0_prec

         ELSEIF ( k .gt. 63 .AND.  myRank == 4 ) THEN

             print*, 'RANK, I,J,K', iTracer, i, j, k
            myFeots % nativeSol % tracer(i,j,k,iTracer) = 1.0_prec

         ENDIF

! Northern boundary

         IF (k .le. 63 .AND. myRank == 5 ) THEN

             print*, 'RANK, I,J,K', iTracer, i, j, k
            myFeots % nativeSol % tracer(i,j,k,iTracer) = 1.0_prec

         ELSEIF ( k .gt. 63 .AND. myRank == 6 ) THEN

             print*, 'RANK, I,J,K', iTracer, i, j, k
            myFeots % nativeSol % tracer(i,j,k,iTracer) = 1.0_prec

         ENDIF


        ENDDO
    !ENDDO
    !ENDIF
 END SUBROUTINE SourceTerms
!  

END PROGRAM FEOTSInitialize
