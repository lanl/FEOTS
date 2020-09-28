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

      
      CALL Init( feots )

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

 SUBROUTINE Init( myFeots )
 ! Sets the source terms and the "relaxation factor" for each tracer
 !
   IMPLICIT NONE
   TYPE( POP_FEOTS ), INTENT(inout) :: myFeots
   ! Local
   INTEGER  :: i, j, k, m, nn, iTracer

      myFeots % nativeSol % tracer = 1.0_prec

 END SUBROUTINE Init
!  

END PROGRAM FEOTSInitialize
