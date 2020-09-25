PROGRAM FEOTSInitialize


! src/POP/
USE POP_FEOTS_Class
USE POP_Native_Class
USE POP_Params_Class


IMPLICIT NONE

   TYPE( POP_FEOTS ) :: feots
   INTEGER :: mpiErr, myRank, nProcs

      CALL MPI_INIT( mpiErr )
      CALL MPI_COMM_RANK( MPI_COMM_WORLD, myRank, mpiErr )
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD, nProcs, mpiErr )

      CALL feots % Build('./runtime.params', myRank, nProcs )

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

      CALL MPI_FINALIZE( mpiErr )

CONTAINS

 SUBROUTINE SourceTerms( myFeots )
 ! Sets the source terms and the "relaxation factor" for each tracer
 !
   IMPLICIT NONE
   TYPE( POP_FEOTS ), INTENT(inout) :: myFeots
   ! Local
   INTEGER  :: i, j, k, m, nn, iTracer

      myFeots % nativeSol % tracer = 1.0_prec

 END SUBROUTINE SourceTerms
!  

END PROGRAM FEOTSInitialize
