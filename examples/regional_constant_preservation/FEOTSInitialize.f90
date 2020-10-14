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

      feots % nativeSol % tracer = 1.0_prec 

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

END PROGRAM FEOTSInitialize
