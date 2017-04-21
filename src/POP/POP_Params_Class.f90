! POP_Params_Class.f90
! 
! Copyright 2016 Joseph Schoonover, Los Alamos National Laboratory (jschoonover@lanl.gov) 
! All rights reserved. 
! 
! POP_Params.f90 is part of the Fast Equilibration of Ocean Tracers Software (FEOTS). 
! 
! //////////////////////////////////////////////////////////////////////////////////////////////// !
 
 
MODULE POP_Params_Class


! src/common/
USE ModelPrecision
USE CommonRoutines
USE ConstantsDictionary
 

 IMPLICIT NONE


    TYPE POP_Params
      ! POPMeshOptions
      INTEGER    :: MeshType 
      INTEGER    :: StencilType
      INTEGER    :: OverlapStencil
      LOGICAL    :: Regional
      REAL(prec) :: south, east, north, west 
      ! TracerModelOptions
      INTEGER    :: TracerModel
      REAL(prec) :: settlingVelocity
      INTEGER    :: nTracers
      INTEGER    :: RunMode
      REAL(prec) :: dt
      INTEGER    :: iterInit
      INTEGER    :: nTimeSteps
      INTEGER    :: nStepsPerDump
      INTEGER    :: nRecordsPerfile
      ! OperatorOptions
      REAL(prec) :: operatorPeriod
      INTEGER    :: nOperatorsPerCycle
      ! FileOptions
      LOGICAL        :: extractRegionalOperators
      CHARACTER(400) :: graphFile
      CHARACTER(400) :: IRFListFile
      INTEGER        :: IRFStart
      CHARACTER(100) :: operatorBaseName
      CHARACTER(400) :: feotsOperatorDirectory
      CHARACTER(400) :: regionalOperatorDirectory
      CHARACTER(400) :: settlingOperatorFile
      INTEGER        :: nIRFFiles
      CHARACTER(400) :: meshfile
      CHARACTER(400) :: regionalMeshfile
      CHARACTER(400) :: outputDirectory
      ! JFNKOptions
      LOGICAL    :: IsPickupRun
      INTEGER    :: maxItersJFNK
      INTEGER    :: maxItersGMRES
      INTEGER    :: mInnerItersGMRES
      INTEGER    :: nResi
      REAL(prec) :: JacobianStepSize
      REAL(prec) :: toleranceJFNK
      REAL(prec) :: toleranceGMRES 
      
       CONTAINS

       PROCEDURE :: Build => Build_POP_Params

    END TYPE POP_Params 
 

 CONTAINS


 SUBROUTINE Build_POP_Params( thisParam )
 ! S/R Build
 !
 !
 ! ==============================================!
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( POP_Params ), intent(out) :: thisParam
   ! LOCAL
   INTEGER :: nUnit
      ! POPMeshOptions
      CHARACTER(20)  :: MeshType 
      CHARACTER(40)  :: StencilType
      LOGICAL        :: Regional
      REAL(prec)     :: south, east, north, west 
      ! TracerModelOptions
      CHARACTER(50)  :: TracerModel
      REAL(prec)     :: settlingVelocity
      INTEGER        :: nTracers
      CHARACTER(20)  :: RunMode
      REAL(prec)     :: dt
      INTEGER        :: iterInit
      INTEGER        :: nTimeSteps
      INTEGER        :: nStepsPerDump
      INTEGER        :: nRecordsPerFile
      ! OperatorOptions
      REAL(prec) :: operatorPeriod
      INTEGER    :: nOperatorsPerCycle
      ! FileOptions
      LOGICAL        :: extractRegionalOperators
      CHARACTER(100) :: operatorBaseName
      CHARACTER(400) :: graphFile
      CHARACTER(400) :: IRFListFile
      INTEGER        :: IRFStart
      CHARACTER(400) :: feotsOperatorDirectory
      CHARACTER(400) :: regionalOperatorDirectory
      CHARACTER(400) :: settlingOperatorFile
      INTEGER        :: nIRFFiles
      CHARACTER(400) :: meshfile
      CHARACTER(400) :: regionalMeshfile
      CHARACTER(400) :: outputDirectory
      ! JFNKOptions
      LOGICAL    :: IsPickupRun
      INTEGER    :: maxItersJFNK
      INTEGER    :: maxItersGMRES
      INTEGER    :: mInnerItersGMRES
      REAL(prec) :: JacobianStepSize
      REAL(prec) :: toleranceJFNK
      REAL(prec) :: toleranceGMRES 

      NAMELIST / POPMeshOptions / MeshType, StencilType, Regional, south, east, north, west
      NAMELIST / TracerModelOptions / TracerModel, settlingVelocity, nTracers, RunMode, dt, iterInit, nTimeSteps, nStepsPerDump, nRecordsPerFile
      NAMELIST / OperatorOptions / operatorPeriod, nOperatorsPerCycle
      NAMELIST / FileOptions / extractRegionalOperators, IRFListFile, IRFStart, feotsOperatorDirectory, &
                               regionalOperatorDirectory, settlingOperatorFile, nIRFFiles, operatorBaseName, &
                               graphFile, regionalMeshFile, meshfile, outputDirectory
      NAMELIST / JFNKOptions / isPickupRun, maxItersJFNK, maxItersGMRES, mInnerItersGMRES, &
                               JacobianStepSize, toleranceJFNK, toleranceGMRES 

      ! Set the default PARAMETERs 
      MeshType         = 'PeriodicTripole'
      StencilType      = 'LaxWendroff'
      Regional         = .TRUE.
      south            = 0.0_prec
      east             = 0.0_prec
      north            = 0.0_prec
      west             = 0.0_prec

      ! TracerModelOptions
      TracerModel      = 'DyeModel'
      settlingVelocity = 0.0_prec
      nTracers         = 1
      RunMode          = 'Forward'
      dt               = 10.0_prec**(-3)
      iterInit         = 0
      nTimeSteps       = 0
      nStepsPerDump    = 0
      nRecordsPerFile  = 10
      ! OperatorOptions
      operatorPeriod     = 86400.0_prec
      nOperatorsPerCycle = 1
      ! FileOptions
      extractRegionalOperators  = .FALSE.
      operatorBaseName          = 'TransportOp'
      graphFile                 = 'graph'
      IRFListFile               = ''
      IRFStart                  = 1
      feotsOperatorDirectory    = './'
      regionalOperatorDirectory = './'
      settlingOperatorFile      = 'settling'
      nIRFFiles                 = 0
      meshfile                  = 'mesh.nc'
      regionalMeshFile          = 'regional_mesh.nc'
      outputdirectory           = './'
      ! JFNKOptions
      isPickupRun      = .FALSE.
      maxItersJFNK     = 500
      maxItersGMRES    = 100
      mInnerItersGMRES = 20
      JacobianStepSize = 10.0_prec**(-2)
      toleranceJFNK    = 10.0_prec**(-7)
      toleranceGMRES   = 10.0_prec**(-7)

      ! Reading in the namelist FILE

      OPEN( UNIT = NEWUNIT(nUnit), FILE = 'runtime.params')
         READ( UNIT = nUnit, NML = POPMeshOptions )
         READ( UNIT = nUnit, NML = TracerModelOptions )
         READ( UNIT = nUnit, NML = OperatorOptions )
         READ( UNIT = nUnit, NML = FileOptions )
         READ( UNIT = nUnit, NML = JFNKOptions )
      CLOSE( UNIT = nUnit ) 

      ! Sanity check - PRINT the results of the namelist READ to screen
     ! WRITE( UNIT = *, NML = POPMeshOptions )
     ! WRITE( UNIT = *, NML = TracerModelOptions )
     ! WRITE( UNIT = *, NML = FileOptions )
     ! WRITE( UNIT = *, NML = JFNKOptions )

      ! Fill in the data structure
      thisParam % MeshType       = GetFlagForChar( TRIM(MeshType) )
      thisParam % StencilType    = GetFlagForChar( TRIM(StencilType) )
      thisParam % Regional       = Regional
      thisParam % south          = south
      thisParam % east           = east
      thisParam % north          = north
      thisParam % west           = west
      ! TracerModelOptions
      thisParam % TracerModel      = GetFlagForChar( TRIM(TracerModel) )
      thisParam % settlingVelocity = settlingVelocity
      !**********************************************************!
      ! If your your model must use a fixed number of tracers,
      ! then you MUST hard-set the "nTracers" attribute here
      ! for that model -- you only need to modify this section
      ! when you are adding a new tracer model that requires
      ! a fixed number of tracers and you want to remove the
      ! control over the choice of the number of tracers from
      ! the scientist's hands.
      IF( thisParam % TracerModel == RadioNuclideModel )THEN 
         thisParam % nTracers = 2
      ELSE
         thisParam % nTracers = nTracers
      ENDIF
      !**********************************************************!

      thisParam % RunMode         = GetFlagForChar( TRIM(RunMode) )
      thisParam % dt              = dt
      thisParam % iterInit        = iterInit
      thisParam % nTimeSteps      = nTimeSteps
      thisParam % nStepsPerDump   = nStepsPerDump
      thisParam % nRecordsPerFile = nRecordsPerFile
      ! OperatorOptions
      thisParam % operatorPeriod     = operatorPeriod
      thisParam % nOperatorsPerCycle = nOperatorsPerCycle
      ! FileOptions
      thisParam % extractRegionalOperators  = extractRegionalOperators
      thisParam % operatorBaseName          = operatorBaseName
      thisParam % graphFile                 = graphFile
      thisParam % feotsOperatorDirectory    = feotsOperatorDirectory
      thisParam % regionalOperatorDirectory = regionalOperatorDirectory
      thisParam % settlingOperatorFile      = settlingOperatorFile
      thisParam % IRFListFile               = IRFListFile
      thisParam % IRFStart                  = IRFStart
      thisParam % nIRFFiles                 = nIRFFiles
      thisParam % outputDirectory           = outputDirectory
      thisParam % regionalMeshFile          = regionalMeshFile
      thisParam % meshfile                  = meshfile
      ! JFNKOptions
      thisParam % isPickupRun      = isPickupRun
      thisParam % maxItersJFNK     = maxItersJFNK
      thisParam % maxItersGMRES    = maxItersGMRES
      thisParam % mInnerItersGMRES = mInnerItersGMRES
      thisParam % nResi            = mInnerItersGMRES*maxItersGMRES
      thisParam % JacobianStepSize = JacobianStepSize
      thisParam % toleranceJFNK    = toleranceJFNK
      thisParam % toleranceGMRES   = toleranceGMRES

 END SUBROUTINE Build_POP_Params

END MODULE POP_Params_Class
