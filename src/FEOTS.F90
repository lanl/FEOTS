PROGRAM FEOTS

USE FEOTS_Driver_Routines

IMPLICIT NONE

  TYPE(FEOTS_CLI) :: cliParams

  CALL cliParams % GetCLIConf( )
  
  IF( cliParams % setupSuccess )THEN

    IF( cliParams % run_Impulse )THEN
      CALL GreedyGraphColoring()
    ELSEIF( cliParams % run_PopMesh )THEN
      CALL GenerateMeshOnlyFile()
    ELSEIF( cliParams % run_genmask )THEN
      CALL GenMask()
    ELSEIF( cliParams % run_regionalExtraction )THEN
      CALL RegionalExtraction(cliParams)
    ELSEIF( cliParams % run_operatorDiagnosis )THEN
      CALL OperatorDiagnosis()
    ELSEIF( cliParams % run_Initialize )THEN
      CALL FEOTSInitialize()
    ELSEIF( cliParams % run_Equilibrator )THEN
      CALL FEOTSEquilibrate()
    ELSEIF( cliParams % run_Integrator )THEN
      CALL FEOTSIntegrate()
    ENDIF

!ExtractOceanState
  ENDIF

END PROGRAM FEOTS
