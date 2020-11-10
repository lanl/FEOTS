! FEOTS_CLI.f90
!
! Author : Joseph Schoonover
! E-mail : jschoonover@lanl.gov, schoonover.numerics@gmail.com
!
! Copyright 2017 Joseph Schoonover, Los Alamos National Laboratory
! 
! Redistribution and use in source and binary forms, with or without
! modification,
! are permitted provided that the following conditions are met:
! 
! 1. Redistributions of source code must retain the above copyright notice, this
! list of conditions and the following disclaimer.
! 
! 2. Redistributions in binary form must reproduce the above copyright notice,
! this list of conditions and the following disclaimer in the 
! documentation and/or other materials provided with the distribution.
! 
! 3. Neither the name of the copyright holder nor the names of its contributors
! may be used to endorse or promote products derived from this 
!  software without specific prior written permission.
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED 
! TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
! PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR 
! CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
! EXEMPLARY,  OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
! BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
! LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE  OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! ////////////////////////////////////////////////////////////////////////////////////////////////

MODULE FEOTS_CLI_Class

USE ModelPrecision
USE CommonRoutines

IMPLICIT NONE

#include "FEOTS_Macros.h"

   TYPE FEOTS_CLI

     LOGICAL :: run_Impulse
     LOGICAL :: run_PopMesh
     LOGICAL :: run_genmask
     LOGICAL :: run_regionalExtraction
     LOGICAL :: run_regionalMaps
     LOGICAL :: run_operatorDiagnosis
     LOGICAL :: run_Initialize
     LOGICAL :: run_Equilibrator
     LOGICAL :: run_Integrator
     LOGICAL :: verticalMixing
     LOGICAL :: setupSuccess
     LOGICAL :: helpNeeded

     LOGICAL :: irfProvided

     CHARACTER(200) :: paramFile
     CHARACTER(200) :: popFile
     CHARACTER(200) :: irfFile
     CHARACTER(200) :: dbRoot
     CHARACTER(200) :: outdir
     CHARACTER(200) :: regionalDb
     INTEGER :: oplevel
     
     CONTAINS

     PROCEDURE :: GetCLIConf
     PROCEDURE :: ValidateCLI
     PROCEDURE :: LogParameters

   END TYPE FEOTS_CLI


CONTAINS

  SUBROUTINE GetCLIConf(cliParams)
    IMPLICIT NONE
    CLASS(FEOTS_CLI), INTENT(out) :: cliParams
    ! Local
    INTEGER :: nArg, argID
    CHARACTER(500) :: argName
    LOGICAL :: helpNeeded, paramFileProvided, irfProvided, oplevelProvided, dbRootProvided, outProvided, mixingProvided, rdbProvided, popFileProvided

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

    ! Default cli parameters
    cliParams % run_Impulse = .FALSE.
    cliParams % run_PopMesh = .FALSE.
    cliParams % run_genmask = .FALSE.
    cliParams % run_regionalExtraction = .FALSE.
    cliParams % run_regionalMaps = .FALSE.
    cliParams % run_operatorDiagnosis = .FALSE.
    cliParams % run_Initialize = .FALSE.
    cliParams % run_Equilibrator = .FALSE.
    cliParams % run_Integrator = .FALSE.
    cliParams % verticalMixing = .TRUE.
    cliParams % helpNeeded = .FALSE.
    cliParams % paramFile = './runtime.params'
    cliParams % irfFile = ''
    cliParams % popFile = ''
    cliParams % oplevel = -1
    cliParams % dbRoot = './'
    cliParams % outdir = './'
    cliParams % regionalDb = './'


    paramFileProvided = .FALSE.
    irfProvided = .FALSE.
    dbRootProvided = .FALSE.
    cliParams % irfProvided = .FALSE.
    oplevelProvided = .FALSE.
    outProvided = .FALSE.


    nArg = command_argument_count( )

    IF( nArg > 0 )THEN
      cliParams % setupSuccess = .TRUE.
      DO argID = 1, nArg

        CALL get_command_argument( argID, argName )

        SELECT CASE( TRIM( argName ) )

          CASE( "impulse" )
            ! GreedyColoring
            cliParams % run_Impulse = .TRUE.

          CASE( "popmesh" )
            !GenerateMeshOnlyFile
            cliParams % run_PopMesh = .TRUE.

          CASE( "genmask" )
            !GenMask
            cliParams % run_genmask = .TRUE.

          CASE( "genmaps" )
            !RegionalMaps 
            cliParams % run_regionalMaps = .TRUE.

          CASE( "region-extraction" )
            !RegionalExtraction 
            cliParams % run_regionalExtraction = .TRUE.

          CASE( "operator-diagnosis" )
            !OperatorDiagnosis
            cliParams % run_operatorDiagnosis = .TRUE.

          CASE( "initialize" )
            !FEOTSInitialize
            cliParams % run_Initialize = .TRUE.

          CASE( "integrate" )
            !FEOTSDriver
            cliParams % run_Integrator = .TRUE.

          CASE( "equilibrate" )
            !FEOTSDriver
            cliParams % run_Equilibrator = .TRUE.

          CASE( "--help" )
            cliParams % helpNeeded = .TRUE.
            EXIT

          CASE( "--param-file" )
            paramFileProvided = .TRUE.

          CASE( "--oplevel" )
            oplevelProvided = .TRUE.

          CASE( "--irf" )
            irfProvided = .TRUE.

          CASE( "--pop-file" )
            popFileProvided = .TRUE.

          CASE( "--dbroot" )
            dbRootProvided = .TRUE.

          CASE( "--out" )
            outProvided = .TRUE.

          CASE( "--regional-db" )
            rdbProvided = .TRUE.

          CASE( "--no-vertical-mixing" )
            cliParams % verticalMixing = .FALSE.


          CASE DEFAULT

            IF( paramFileProvided )THEN

              cliParams % paramFile = TRIM( argName )
              paramFileProvided = .FALSE.

            ENDIF

            IF( oplevelProvided )THEN

              READ(argName,*)  cliParams % opLevel
              opLevelProvided = .FALSE.

            ENDIF

            IF( irfProvided )THEN

              cliParams % irfFile = TRIM( argName )
              irfProvided = .FALSE.
              cliParams % irfProvided = .TRUE.

            ENDIF

            IF( popFileProvided )THEN

              cliParams % popFile = TRIM( argName )
              popFileProvided = .FALSE.

            ENDIF

            IF( dbRootProvided )THEN

              cliParams % dbRoot = TRIM( argName )
              dbRootProvided = .FALSE.

            ENDIF

            IF( outProvided )THEN

              cliParams % outdir = TRIM( argName )
              outProvided = .FALSE.

            ENDIF

            IF( rdbProvided )THEN

              cliParams % regionalDb = TRIM( argName )
              rdbProvided = .FALSE.

            ENDIF

        END SELECT

      ENDDO

    ELSE

      cliParams % helpNeeded = .TRUE.

    ENDIF

    IF( cliParams % helpNeeded ) THEN

      cliParams % setupSuccess = .FALSE.

      PRINT*, '  feots [tool] [options]'      
      PRINT*, ' '
      PRINT*, ' [tool] can be :'
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
      PRINT*, '     You must specify the IRF file using the --irf option.'
      PRINT*, ' '
      PRINT*, '   region-extraction'
      PRINT*, '     Create regional transport operators from global transport operators. Regional'
      PRINT*, '     operators are stored in the --regional-db directory.'
      PRINT*, ' '
      PRINT*, '   genmaps'
      PRINT*, '     Create a mappings.regional file from a valid mask file. The mappings.regional'
      PRINT*, '     file is stored in the --out directory.'
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
      PRINT*, '   --help'
      PRINT*, '     Display this help message'
      PRINT*, ' '
      PRINT*, '    --param-file /path/to/param/file'
      PRINT*, '       Specifies the full path to a file with namelist settings for'
      PRINT*, '       the feots application. If not provided, runtime.params in  '
      PRINT*, '       your current directory is assumed.                          '
      PRINT*, ' '
      PRINT*, '    --pop-file /path/to/irf-file'
      PRINT*, '       Specifies the full path to a netcdf file with standard POP output'
      PRINT*, '       (For popmesh)'
      PRINT*, ' '
      PRINT*, '    --irf /path/to/irf-file'
      PRINT*, '       Specifies the full path to a netcdf file with IRFs'
      PRINT*, '       (For operator diagnosis and regional extraction)'
      PRINT*, ' '
      PRINT*, '    --oplevel 0'
      PRINT*, '       Specifies the index of the operator in the operator sequence'
      PRINT*, '       This option determines the time level encoded to _advect.{oplevel}.data/conn'
      PRINT*, ' '
      PRINT*, '    --dbroot /path/to/feot/db'
      PRINT*, '       Specifies the path to a FEOTS database'
      PRINT*, ' '
      PRINT*, '    --out /path/to/output/directory'
      PRINT*, '       Specifies the path to write model output. Defaults to ./'
      PRINT*, ' '
      PRINT*, '    --no-vertical-mixing'
      PRINT*, '       Disables the vertical mixing operator for forward integration and equilibration'
      PRINT*, ' '
      PRINT*, '    --regional-db /path/to/regional-database/directory'
      PRINT*, '       Specifies the path to read/write regional operators. Defaults to ./'
      PRINT*, ' '
      PRINT*, ' '
      PRINT*, ' '

      RETURN

    ENDIF

    CALL cliParams % ValidateCLI( )

    CALL cliParams % LogParameters( )

  END SUBROUTINE GetCLIConf

  SUBROUTINE LogParameters(cliParams)
#undef __FUNC__
#define __FUNC__ "LogParameters"
    IMPLICIT NONE
    CLASS(FEOTS_CLI), INTENT(in) :: cliParams

      INFO('paramFile = '//TRIM(cliParams % paramFile))
      INFO('irfFile = '//TRIM(cliParams % irfFile))
      INFO('dbRoot = '//TRIM(cliParams % dbRoot))
      INFO('outdir = '//TRIM(cliParams % outdir))
      INFO('regionalDb = '//TRIM(cliParams % regionalDb))

  END SUBROUTINE LogParameters

  SUBROUTINE ValidateCLI(cliParams)
    IMPLICIT NONE
    CLASS(FEOTS_CLI), INTENT(inout) :: cliParams


      cliParams % setupSuccess = .TRUE.
      IF( cliParams % helpNeeded )THEN
        cliParams % setupSuccess = .FALSE.
      ELSE
    
        IF( cliParams % run_operatorDiagnosis )THEN
          IF( TRIM(cliParams % irfFile) == '' )THEN
            PRINT*, 'ERROR: IRF File needed for operator-diagnosis.'
            cliParams %  setupSuccess = .FALSE.
          ENDIF

          IF( cliParams % opLevel < 0 )THEN
            PRINT*, 'WARNING: Operator level not provided. Transport operators will use sequence level of 0'
            cliParams % opLevel = 0
            cliParams %  setupSuccess = .TRUE.
          ENDIF
        ENDIF
        IF( cliParams % run_regionalExtraction )THEN

          IF( cliParams % opLevel < 0 )THEN
            PRINT*, 'WARNING: Operator level not provided. Transport operators will use sequence level of 0'
            cliParams % opLevel = 0
            cliParams %  setupSuccess = .TRUE.
          ENDIF
        ENDIF

      ENDIF

  END SUBROUTINE ValidateCLI
        

END MODULE FEOTS_CLI_Class
