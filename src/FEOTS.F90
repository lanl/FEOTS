PROGRAM FEOTS

USE FEOTS_Driver_Routines

IMPLICIT NONE

  LOGICAL :: run_Impulse
  LOGICAL :: run_PopMesh
  LOGICAL :: run_genmask
  LOGICAL :: run_regionalExtraction
  LOGICAL :: run_operatorDiagnosis
  LOGICAL :: run_Initialize
  LOGICAL :: run_Equilibrator
  LOGICAL :: run_Integrator
  LOGICAL :: setupSuccess
  CHARACTER(500) :: paramFile

  CALL GetCLIConf( )
  
  IF( setupSuccess )THEN

    IF( run_Impulse )THEN
      CALL GreedyGraphColoring()
    ELSEIF( run_PopMesh )THEN
      CALL GenerateMeshOnlyFile()
    ELSEIF( run_genmask )THEN
      CALL GenMask()
    ELSEIF( run_regionalExtraction )THEN
      CALL RegionalExtraction()
    ELSEIF( run_operatorDiagnosis )THEN
      CALL OperatorDiagnosis()
    ELSEIF( run_Initialize )THEN
      CALL FEOTSInitialize()
    ELSEIF( run_Equilibrator )THEN
      CALL FEOTSEquilibrate()
    ELSEIF( run_Integrator )THEN
      CALL FEOTSIntegrate()
    ENDIF

!ExtractOceanState
  ENDIF

CONTAINS

  SUBROUTINE GetCLIConf( )
    
    ! Local
    INTEGER :: nArg, argID
    CHARACTER(500) :: argName
    LOGICAL :: helpNeeded, paramFileProvided

    run_Impulse = .FALSE.
    run_PopMesh = .FALSE.
    run_genmask = .FALSE.
    run_regionalExtraction = .FALSE.
    run_operatorDiagnosis = .FALSE.
    run_Initialize = .FALSE.
    run_Equilibrator = .FALSE.
    run_Integrator = .FALSE.
    helpNeeded = .FALSE.
    paramFileProvided = .FALSE.

    paramFile = './runtime.params'

    nArg = command_argument_count( )

    IF( nArg > 0 )THEN
      setupSuccess = .TRUE.
      DO argID = 1, nArg

        CALL get_command_argument( argID, argName )

        SELECT CASE( TRIM( argName ) )

          CASE( "impulse" )
            ! GreedyColoring
            run_Impulse = .TRUE.
            setupSuccess = .TRUE.

          CASE( "popmesh" )
            !GenerateMeshOnlyFile
            run_PopMesh = .TRUE.
            setupSuccess = .TRUE.

          CASE( "genmask" )
            !GenMask
            run_genmask = .TRUE.
            setupSuccess = .TRUE.

          CASE( "region-extraction" )
            !RegionalExtraction 
            run_regionalExtraction = .TRUE.
            setupSuccess = .TRUE.

          CASE( "operator-diagnosis" )
            !OperatorDiagnosis
            run_operatorDiagnosis = .TRUE.
            setupSuccess = .TRUE.

          CASE( "initialize" )
            !FEOTSInitialize
            run_Initialize = .TRUE.
            setupSuccess = .TRUE.

          CASE( "integrate" )
            !FEOTSDriver
            run_Integrator = .TRUE.
            setupSuccess = .TRUE.

          CASE( "equilibrate" )
            !FEOTSDriver
            run_Equilibrator = .TRUE.
            setupSuccess = .TRUE.

          CASE( "help" )
            helpNeeded   = .TRUE.
            setupSuccess = .FALSE.

          CASE( "--param-file" )
            paramFileProvided = .TRUE.

          CASE DEFAULT

            setupSuccess      = .TRUE.

            IF( paramFileProvided )THEN

              paramFile = TRIM( argName )
              paramFileProvided = .FALSE.

            ENDIF

        END SELECT

      ENDDO

    ELSE

      helpNeeded = .TRUE.
      setupSuccess = .FALSE.

    ENDIF

    IF( helpNeeded ) THEN

      PRINT*, 'FEOTS (feots) Command Line Interface' 
      PRINT*, ' Copyright Los Alamos National Laboratory (2017-2020)'
      PRINT*, ' Licensed for use under 3-Clause BSD License'
      PRINT*, ' '
      PRINT*, ' For support related issues, https://github.com/lanl/feots/issues/new'
      PRINT*, ' '
      PRINT*, ' A program for performing creating impulse functions, diagnosing transport'
      PRINT*, ' operators from POP IRFs, and conducting offline tracer simulations using '
      PRINT*, ' diagnosed transport operators.'
      PRINT*, ' '
      PRINT*, '  feots [tool] [options]'      
      PRINT*, ' '
      PRINT*, ' [tool] can be :'
      PRINT*, ' '
      PRINT*, '   help'
      PRINT*, '     Display this help message'
      PRINT*, ' '
      PRINT*, '   impulse'
      PRINT*, '     Use a POP-Mesh, with land-mask, and a chosen advection-difussion stencil'
      PRINT*, '     to create impulse fields for capturing impulse response functions.'
      PRINT*, ' '
      PRINT*, '   popmesh'
      PRINT*, '     Extract POP-Mesh information from POP standard output.'
      PRINT*, ' '
      PRINT*, '   genmask'
      PRINT*, '     Create a regional FEOTS mask using lat-lon domain bounds'
      PRINT*, ' '
      PRINT*, '   operator-diagnosis'
      PRINT*, '     Diagnose transport operators using impulse fields and POP IRF output.'
      PRINT*, ' '
      PRINT*, '   region-extraction'
      PRINT*, '     Create regional transport operators from global transport operators'
      PRINT*, ' '
      PRINT*, '   initialize'
      PRINT*, '     Use the built in initialization routines to create tracer initial conditions'
      PRINT*, ' '
      PRINT*, '   integrate'
      PRINT*, '     Run the offline tracer simulation in a forward integration mode'
      PRINT*, ' '
      PRINT*, '   equilibrate'
      PRINT*, '     Run the offline tracer simulation using JFNK to find the equilibrated tracer field'
      PRINT*, ' '
      PRINT*, '  [options] can be :'
      PRINT*, ' '
      PRINT*, '    --param-file /path/to/param/file'
      PRINT*, '       Specifies the full path to a file with namelist settings for'
      PRINT*, '       the feots application. If not provided, runtime.params in  '
      PRINT*, '       your current directory is assumed.                          '
      PRINT*, ' '
      PRINT*, ' '

      setupSuccess = .FALSE.
      RETURN
    ENDIF



  END SUBROUTINE GetCLIConf

END PROGRAM FEOTS
