PROGRAM FEOTS

USE FEOTS_Driver_Routines

IMPLICIT NONE

  TYPE(FEOTS_CLI) :: cliParams

  CALL cliParams % GetCLIConf( )

#ifdef _OPENMP
!$OMP PARALLEL
   IF(OMP_GET_THREAD_NUM() == 0)THEN
     PRINT*, "Number of OpenMP threads = ", OMP_GET_NUM_THREADS()
   ENDIF
!$OMP END PARALLEL
#endif
  
  IF( cliParams % setupSuccess )THEN

    IF( cliParams % run_Impulse )THEN

      CALL GreedyGraphColoring(cliParams)

    ELSEIF( cliParams % run_PopMesh )THEN

      CALL GenerateMeshOnlyFile(cliParams)

    ELSEIF( cliParams % run_genmask )THEN

      CALL GenMask(cliParams)

    ELSEIF( cliParams % run_regionalExtraction )THEN

      CALL RegionalExtraction(cliParams)

    ELSEIF( cliParams % run_operatorDiagnosis )THEN

      CALL OperatorDiagnosis(cliParams)

    ELSEIF( cliParams % run_Initialize )THEN

      CALL FEOTSInitialize(cliParams)

    ELSEIF( cliParams % run_Equilibrator )THEN

      CALL FEOTSEquilibrate(cliParams)

    ELSEIF( cliParams % run_Integrator )THEN

      CALL FEOTSIntegrate(cliParams)

    ENDIF

  ENDIF

END PROGRAM FEOTS
