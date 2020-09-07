! POP_FEOTS_Class.f90
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
! AND
! ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED 
! TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
! PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR 
! CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
! EXEMPLARY,
! OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
! BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
! LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE
! OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !
 
 
 MODULE POP_FEOTS_Class

! src/common/
USE ModelPrecision
USE ConstantsDictionary
! src/matrices
USE CRSMatrix_Class
! src/solutionstorage/
USE TracerStorage_Class
! src/POP/
USE POP_Params_Class
USE POP_Stencil_Class
USE POP_Mesh_Class
USE POP_Regional_Class
USE POP_Native_Class
USE FEOTS_CLI_Class

 IMPLICIT NONE


! ==================================== POP_FEOTS Description ====================================== !
!
! 18.04.2016
!
! The POP_FEOTS Class was written in order to facilitate explicit integration of the TracerStorage
! Class and provide routines for file I/O in binary and tecplot formats. This class combines the
! POPMesh and TracerStorage classes to make it easy to setup and run a transport problem
! with diagnosed transport operators.
!
! Support routines are provided to simplify the setup of a problem. Relaxation source terms and 
! "hard-set" fields can be constructed using 3-D arrays that correspond to the POP mesh storage
! structure. The 3-D fields can be passed to a routine "MapToDOF" to convert the 3-D array into
! the 1-D array structure that is used to conduct the forward integration and/or iterative solve.
! Additionally "MapFromDOF" maps the 1-D array storage back to the POP 3-D array to prepare for
! file output ( to TecPlot or binary ).
! ================================================================================================ !
!

! Data Storage requirements
!  solution 
!  > 5*INTEGER ( nDOF, nPeriods, nOps, nTracers, currentPeriod )
!  > 2*prec    ( opPeriod, dt )
!  > CRSMatrix ( transportOp )
!  > > 3*INTEGER (nRows, nCols, nElems)
!  > > 9*nDOF*prec ( A )
!  > > 2*nDOF*INTEGER ( rowBounds )
!  > > nElems*INTEGER ( col )
!
!  > CRSMatrix ( diffusionOp )
!  > > 3*INTEGER (nRows, nCols, nElems)
!  > > 3*nDOF*prec ( A )
!  > > 2*nDOF*INTEGER ( rowBounds )
!  > > nElems*INTEGER ( col )
!
!  > 4*nDOF*nTracers*prec ( tracers, source, rFac, mask )
!  > nDOF*prec ( volume )
!  > nTracers*INTEGER (tracerIDs)
!
!  feotsMap
!  > 4*prec (south, north, east, west)
!  > LOGICAL (crossesPrimeMeridian)
!  > 3*INTEGER (nMasks, nCells, nDOF)
!  > 3*nDOF (ijkInRegion)
!  > nDOF (ijkInRegion)
!  > nGlobalDOF (inverseDOFMap)
!  > 3*nDOF (dofToLocalIJK)
!  > nTracers*BoundaryMap
!  > > 2*INTEGER (nBCells, nPCells)
!  > > nBCells (boundaryCells)
!  > > nPCells (prescribedCells)
!
!
!  mesh
!  5*INTEGER ( nX,nY,nZ,nDOF,meshType )
!  nX*nY*prec (tLon)
!  nX*nY*prec (tLat)
!  nX*nY*prec (dXt)
!  nX*nY*prec (dYt)
!  nX*nY*prec (tArea)
!  nX*nY*prec (KmT)
!  nZ*prec (z)
!  nZ*prec (dz)
!  nZ*prec (dzw)
!  nX*nY*nZ*prec (tracermask)
!  3*nDOF*INTEGER (DOFtoIJK)
!  nX*nY*nZ*INTEGER (IJKtoDOF)
!
!  nativeSol
!  4*INTEGER ( nX, nY, nZ, nTracers )
!  nTracers*INTEGER ( tracerIds )
!  nX*nY*nZ*nTracers*prec ( tracer )
!  nX*nY*nZ*nTracers*prec ( mask )
!  nX*nY*nZ*nTracers*prec ( source )
!  nX*nY*nZ*nTracers*prec ( rfac )
!  nX*nY*nZ*prec ( volume )
!  nX*nY*nZ*prec ( temperature )
!  nX*nY*nZ*prec ( salinity )
!  nX*nY*nZ*prec ( density )

   TYPE POP_FEOTS
     
      TYPE( TracerStorage ) :: solution
      TYPE( POP_Regional ) :: feotsMap
      TYPE( POP_Mesh ) :: mesh
      TYPE( POP_Params ) :: params
      TYPE( POP_Native ) :: nativeSol
      INTEGER :: myRank, nProcs

      !! Water Mass Tagging
      !REAL(prec)                :: stateMask(1:3)
      !REAL(prec), ALLOCATABLE   :: stateLowerBound(:)
      !REAL(prec), ALLOCATABLE   :: stateUpperBound(:)

      CONTAINS

      PROCEDURE :: Build => Build_POP_FEOTS
      PROCEDURE :: Trash => Trash_POP_FEOTS

      PROCEDURE :: MapAllToDOF      => MapAllToDOF_POP_FEOTS
      PROCEDURE :: MapTracerToDOF   => MapTracerToDOF_POP_FEOTS
      PROCEDURE :: MapHSetToDOF     => MapHSetToDOF_POP_FEOTS
      PROCEDURE :: MapSourceToDOF   => MapSourceToDOF_POP_FEOTS
      PROCEDURE :: MapTracerFromDOF => MapTracerFromDOF_POP_FEOTS

      PROCEDURE :: LoadNewStates       => LoadNewStates_POP_FEOTS
     
      PROCEDURE :: StepForward
      PROCEDURE :: ForwardStep         => ForwardStep_POP_FEOTS
      PROCEDURE :: ForwardStepEuler    => ForwardStepEuler_POP_FEOTS
      PROCEDURE :: ForwardStepAB2      => ForwardStepAB2_POP_FEOTS
      PROCEDURE :: ForwardStepAB3      => ForwardStepAB3_POP_FEOTS

      PROCEDURE :: CycleIntegration      => CycleIntegration_POP_FEOTS
      PROCEDURE :: CycleIntegrationEuler => CycleIntegrationEuler_POP_FEOTS
      PROCEDURE :: CycleIntegrationAB2   => CycleIntegrationAB2_POP_FEOTS
      PROCEDURE :: CycleIntegrationAB3   => CycleIntegrationAB3_POP_FEOTS
      PROCEDURE :: DotProduct            => DotProduct_POP_FEOTS
      PROCEDURE :: JacobianAction        => JacobianAction_POP_FEOTS
      PROCEDURE :: SolveGMRES            => SolveGMRES_POP_FEOTS
      PROCEDURE :: JFNK                  => JFNK_POP_FEOTS

      PROCEDURE :: VerticalMixing => VerticalMixing_POP_FEOTS
      PROCEDURE :: VerticalMixingAction
      PROCEDURE :: GetMax
      PROCEDURE :: Mask

   END TYPE POP_FEOTS

   REAL(prec), PARAMETER, PRIVATE :: cg_tolerance = 1.0D-7
   INTEGER, PARAMETER, PRIVATE    :: cg_itermax   = 50000

CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE Build_POP_FEOTS( this, cliParams, myRank, nProcs ) 
 ! S/R Build
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( POP_FEOTS ), INTENT(out) :: this
   TYPE( FEOTS_CLI ) :: cliParams
   INTEGER, INTENT(in)             :: myRank, nProcs
   ! Local
   INTEGER    :: nX, nY, nZ, nTracers, nRow, nPeriods
   INTEGER    :: nDOF, iLayer, iMask, iTracer
   REAL(prec) :: opPeriod, dt
   INTEGER    :: nOps, i, j, k, m, stencilSize, fUnit
   INTEGER    :: trackingVar, tracerID
   CHARACTER(200) :: oceanStateFile
   CHARACTER(20)  :: trackingVar_Char
   CHARACTER(5)   :: fileIDChar
   

      PRINT*, 'S/R : Build_POP_FEOTS : Start...'
      this % myRank = myRank
      this % nProcs = nProcs

      CALL this % params % Build(TRIM(cliParams % paramFile))
      this % params % dbRoot = cliParams % dbRoot

      IF( this % params % TracerModel /= DyeModel )THEN
         PRINT*, 'MPI currently only configured for Passive Dye Model.'
         STOP 'Stopping!'
      ENDIF

      IF( nProcs /= this % params % nTracers )THEN
         PRINT*, 'Number of MPI ranks must equal the number of tracers.', this % params % nTracers
         STOP 'Stopping!'
      ENDIF


      IF( this % params % StencilType == LaxWendroff )THEN
         stencilSize = 9
      ELSEIF( this % params % StencilType == LaxWendroff27 )THEN
         stencilSize = 27
      ELSEIF( this % params % StencilType == Upwind3 )THEN
         stencilSize = 13
      ELSE
         STOP 'Bad Stencil!'
      ENDIF

      IF( this % params % Regional )THEN
         CALL this % mesh % LoadWithMask( TRIM(this % params % RegionalMeshFile) )
         CALL this % feotsMap % ReadPickup( TRIM(this % params % regionalOperatorDirectory)//'mappings', maskProvided=.TRUE. )

         !!this % mesh % DOFtoIJK = this % feotsMap % dofToLocalIJK
         DO m = 1, this % feotsMap % nCells
            i = this % feotsMap % dofToLocalIJK(1,m)
            j = this % feotsMap % dofToLocalIJK(2,m)
            k = this % feotsMap % dofToLocalIJK(3,m)
            this % mesh % IJKtoDOF(i,j,k) = m
            this % mesh % DOFtoIJK(1,m)   = i
            this % mesh % DOFtoIJK(2,m)   = j
            this % mesh % DOFtoIJK(3,m)   = k
         ENDDO
         nDOF = this % feotsMap % nCells
         this % mesh % nDOF = nDOF
      ELSE
         CALL this % mesh % Load( TRIM(this % params % dbRoot)//'/mesh/mesh.nc' )
         nDOF = this % mesh % nDOF
      ENDIF

      ! Allocates space for the solution storage as a 1-D array and allocates
      ! space for the transport operators 
      CALL this % solution % Build( nDOF, nOps, (/nDOF*stencilSize,nDOF*3/), &
                                    this % params % nOperatorsPerCycle, &
                                    this % params % nTracers, & 
                                    this % params % operatorPeriod, &
                                    this % params % dt, &
                                    myRank, nProcs )

      CALL this % nativeSol % Build( this % mesh, this % params % nTracers, myRank, nProcs )

      IF( this % params % Regional )THEN
        DO iTracer = 1, this % nativeSol  % nTracers
          tracerID = this % nativeSol % tracerIds(iTracer)
          DO m = 1, this % feotsMap % bMap(tracerID) % nBCells
            i = this % feotsMap % dofToLocalIJK(1,this % feotsMap % bMap(tracerID) % boundaryCells(m))
            j = this % feotsMap % dofToLocalIJK(2,this % feotsMap % bMap(tracerID) % boundaryCells(m))
            k = this % feotsMap % dofToLocalIJK(3,this % feotsMap % bMap(tracerID) % boundaryCells(m))
            this % nativeSol % mask(i,j,k,iTracer) = 0.0_prec
          ENDDO
        ENDDO
      ENDIF
!

!      IF( this % params % WaterMassTagging )THEN
!         
!         PRINT*, '   Enabling water mass tagging.'
!         ! Overwrite the number of tracers
!         this % params % nTracers = this % params % nLayers*this % feotsMap % nMasks
!
!         ALLOCATE( this % stateLowerBound(1:this % params % nLayers), &
!                   this % stateUpperBound(1:this % params % nLayers) )
!         this % stateLowerBound = 0.0_prec
!         this % stateUpperBound = 0.0_prec
!         
!         !IF( LoadForDriver )THEN
!            OPEN( UNIT   = NewUnit( fUnit ), &
!                  FILE   = 'watermass.config', &
!                  FORM   = 'FORMATTED', &
!                  STATUS = 'OLD', &
!                  ACTION = 'READ' )
!
!            READ( fUnit, '(A20)' ) trackingVar_Char
!            PRINT*, '   Tagging water masses with '//TRIM( trackingVar_Char )
!            trackingVar = GetFlagforChar( TRIM(trackingVar_Char) ) 
!            IF( trackingVar == Temperature ) THEN
!               this % stateMask(1) = 1.0_prec
!               this % stateMask(2) = 0.0_prec
!               this % stateMask(3) = 0.0_prec
!            ELSEIF( trackingVar == Salinity ) THEN
!               this % stateMask(1) = 0.0_prec
!               this % stateMask(2) = 1.0_prec
!               this % stateMask(3) = 0.0_prec
!            ELSEIF( trackingVar == Density ) THEN
!               this % stateMask(1) = 0.0_prec
!               this % stateMask(2) = 0.0_prec
!               this % stateMask(3) = 1.0_prec
!            ELSE
!               PRINT*, 'Bad Water Mass Tracking variable.'
!               STOP 'STOPPING!'
!            ENDIF
!
!            PRINT*, '         Layer |        Lower Bound           |    Upper Bound' 
!            DO i = 1, this % params % nLayers
!               READ( fUnit, * ) this % stateLowerBound(i), this % stateUpperBound(i)
!               PRINT*, i, '   |', this % stateLowerBound(i),'   |', this % stateUpperBound(i)
!            ENDDO
!
!            CLOSE( fUnit )
!         !ENDIF
!      ENDIF
!
!      IF( this % params % TracerModel == RadionuclideModel .OR. &
!          this % params % TracerModel == SettlingModel )THEN 
!         PRINT*, '  Setting up Settling Operator'
!         CALL this % SetupSettlingOperator( )
!      ENDIF

      ! Allocates space for the solution storage on the native mesh 
!#ifdef HAVE_MPI
!      CALL this % nativeSol % Build( this % mesh, 1, myRank )
!#else
      !ELSE
  

  ! IF( this % params % waterMassTagging )THEN

  !    WRITE(fileIDChar, '(I5.5)' ) 1
  !    oceanStateFile = TRIM(this % params % regionalOperatorDirectory)//'Ocean.'//fileIDChar//'.nc'
  !    PRINT*, 'Loading initial ocean state : ', TRIM(oceanStateFile)
  !    CALL this % nativeSol % LoadOceanState( this % mesh, TRIM(oceanStateFile) )

  !    DO iMask = 1, this % feotsMap % nMasks
  !       DO iLayer = 1, this % params % nLayers
  !             
  !          iTracer = iLayer + (iMask-1)*( this % params % nLayers )
  !          DO m = 1, this % feotsMap % bMap(iMask) % nBCells
  !             i = this % feotsMap % dofToLocalIJK(1,this % feotsMap % bMap(iMask) % boundaryCells(m))
  !             j = this % feotsMap % dofToLocalIJK(2,this % feotsMap % bMap(iMask) % boundaryCells(m))
  !             k = this % feotsMap % dofToLocalIJK(3,this % feotsMap % bMap(iMask) % boundaryCells(m))
  !             this % nativeSol % mask(i,j,k,iTracer) = 0.0_prec
  !          ENDDO
  !       ENDDO
  !    ENDDO

  ! ENDIF

#ifdef HAVE_MPI
      ! Reset the number of tracers (for each rank) to 1
      this % params % nTracers = 1
#endif
      PRINT*, 'S/R : Build_POP_FEOTS : Finish.'

 END SUBROUTINE Build_POP_FEOTS
!
 SUBROUTINE Trash_POP_FEOTS( this )
 ! S/R Trash
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( POP_FEOTS ), INTENT(inout) :: this

      CALL this % mesh % Trash( )
      CALL this % feotsMap % Trash( )
      CALL this % solution % Trash( )
      CALL this % nativeSol % Trash( )     

 END SUBROUTINE Trash_POP_FEOTS
!
! SUBROUTINE SetupSettlingOperator_POP_FEOTS( this )
!   IMPLICIT NONE
!   CLASS( POP_FEOTS ), INTENT(inout) :: this
!   ! Local
!   REAL(prec) :: ws, fac
!   INTEGER    :: iel, row, col, i, j, k
!
!      ws = this % params % settlingVelocity
!
!      iel = 0
!      DO row = 1, this % mesh % nDOF
!
!         i = this % mesh % DOFtoIJK(1,row) ! 
!         j = this % mesh % DOFtoIJK(2,row) ! 
!         k = this % mesh % DOFtoIJK(3,row) ! vertical level
!         
!         IF( k > 1 )THEN
!
!            fac = ws/( this % mesh % z(k) - this % mesh % z(k-1) )
!          
!            col = this % mesh % IJKtoDOF(i,j,k-1)
!            iel = iel + 1
!            this % solution % transportOps(3) % rowBounds(1,row) = iel
!            this % solution % transportOps(3) % A(iel)           = fac              
!            this % solution % transportOps(3) % col(iel)         = col ! sub-diagonal   
! 
!            iel = iel + 1
!            this % solution % transportOps(3) % rowBounds(2,row) = iel
!            this % solution % transportOps(3) % A(iel)           = -fac              
!            this % solution % transportOps(3) % col(iel)         = row ! diagonal              
! 
!         ELSE
!
!            ! Here, k=1, ie we're at the top-most layer. The flux through the
!            ! face is zero
!            fac = ws/( this % mesh % z(k)  )
!            iel = iel + 1
!            this % solution % transportOps(3) % rowBounds(1,row) = iel
!            this % solution % transportOps(3) % rowBounds(2,row) = iel
!            this % solution % transportOps(3) % A(iel)           = -fac              
!            this % solution % transportOps(3) % col(iel)         = row ! diagonal              
! 
!         ENDIF
!
!      ENDDO
!
! END SUBROUTINE SetupSettlingOperator_POP_FEOTS
!
!
!==================================================================================================!
!----------------------------------------- Type Specific ------------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE MapAllToDOF_POP_FEOTS( this )
 ! S/R MapAllToDOF
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( POP_FEOTS ), INTENT(inout) :: this
   
      CALL this % MapTracerToDOF( )
      CALL this % MapHSetToDOF( )
      CALL this % MapSourceToDOF( )

 END SUBROUTINE MapAllToDOF_POP_FEOTS
!
 SUBROUTINE MapTracerToDOF_POP_FEOTS( this )
 ! S/R MapTracerToDOF
 !
 !   This subroutine maps the attribute "nativeTracer" to the DOF array in the attribute 
 !   " solution % tracers ". It is assumed that "nativeTracer" is filled in with the desired tracer
 !   field before calling this routine. 
 !
 ! =============================================================================================== !
   IMPLICIT NONE 
   CLASS( POP_FEOTS ), INTENT(inout) :: this
   ! Local
   INTEGER :: i

      DO i = 1, this % solution % nTracers 
         this % solution % tracers(:,i) = this % mesh % MapFromIJKtoDOF( this % nativeSol % tracer(:,:,:,i) )
      ENDDO
      this % solution % volume = this % mesh % MapFromIJKtoDOF( this % nativeSol % volume )
   
 END SUBROUTINE MapTracerToDOF_POP_FEOTS
!
 SUBROUTINE MapHSetToDOF_POP_FEOTS( this )
 ! S/R MapHSetToDOF
 !
 !   This subroutine maps the attribute "nativeHSet" and "nativeMask" to the DOF array in the 
 !   attributes "solution % hSetTracers" and "solution % mask". It is assumed that "nativeHSet" and 
 !   "nativeMask" are filled in with the desired hard-set tracer field and hard-set mask before 
 !   calling this routine. 
 !
 ! =============================================================================================== !
   IMPLICIT NONE 
   CLASS( POP_FEOTS ), INTENT(inout) :: this
   ! Local
   INTEGER :: i

      DO i = 1, this % solution % nTracers
         this % solution % mask(:,i) = this % mesh % MapFromIJKtoDOF( this % nativeSol % mask(:,:,:,i) )
      ENDDO
   
 END SUBROUTINE MapHSetToDOF_POP_FEOTS
!
 SUBROUTINE MapSourceToDOF_POP_FEOTS( this )
 ! S/R MapSourceToDOF
 !
 !   This subroutine maps the attribute "nativeSource" and "nativeRFac" to the DOF array in the 
 !   attributes "solution % source" and "solution % rFac". It is assumed that "nativeSource" and 
 !   "nativeRFac" are filled in with the desired source field and relaxation factor before 
 !   calling this routine. 
 !
 ! =============================================================================================== !
   IMPLICIT NONE 
   CLASS( POP_FEOTS ), INTENT(inout) :: this
   ! Local
   INTEGER :: i

      DO i = 1, this % solution % nTracers
         this % solution % source(:,i) = this % mesh % MapFromIJKtoDOF( this % nativeSol % source(:,:,:,i) )
         this % solution % rFac(:,i) = this % mesh % MapFromIJKtoDOF( this % nativeSol % rFac(:,:,:,i) )
      ENDDO
   
 END SUBROUTINE MapSourceToDOF_POP_FEOTS
!
 SUBROUTINE MapTracerFromDOF_POP_FEOTS( this )
 ! S/R MapTracerFromDOF
 !
 !   This subroutine maps the attribute "solution % tracers" to the POP-native array in the attribute 
 !   "nativeTracer ".
 !
 ! =============================================================================================== !
   IMPLICIT NONE 
   CLASS( POP_FEOTS ), INTENT(inout) :: this
   ! Local
   INTEGER :: i

      DO i = 1, this % solution % nTracers
         this % nativeSol % tracer(:,:,:,i) = this % mesh % MapFromDOFtoIJK( this % solution % tracers(:,i) )
      ENDDO
      this % nativeSol % volume = this % mesh % MapFromDOFtoIJK( this % solution % volume )
   
 END SUBROUTINE MapTracerFromDOF_POP_FEOTS
!
 SUBROUTINE LoadNewStates_POP_FEOTS( this, tn, myRank, operatorPeriod )
   IMPLICIT NONE
   CLASS( POP_FEOTS ), INTENT(inout) :: this
   REAL(prec), INTENT(in)            :: tn
   INTEGER, INTENT(in)               :: myRank
   INTEGER, INTENT(in), OPTIONAL     :: operatorPeriod
   ! Local
   INTEGER         :: i, j, k, m, iT, dof, iTracer, iLayer, iMask
   INTEGER :: row, col, ii, jj, kk
   LOGICAL         :: operatorsSwapped
   CHARACTER(400)  :: oceanStateFile
   CHARACTER(500)  :: fileBase
   CHARACTER(5)    :: fileIDChar
   REAL(prec)      :: trackingVar
   
      !$OMP MASTER
      IF( PRESENT(operatorPeriod) )THEN
         
         WRITE( fileIDChar, '(I5.5)' ) operatorPeriod
         IF( this % params % Regional )THEN
            fileBase = TRIM(this % params % regionalOperatorDirectory)
         ELSE
            fileBase = TRIM(this % params % dbRoot)//'/ops'
         ENDIF
         PRINT*, '  Loading Operator : '//TRIM(fileBase)//'/transport.'//fileIDChar//'.h5'
         CALL this % solution % transportOp % ReadCRSMatrix_HDF5( TRIM(fileBase)//'/transport.'//fileIDChar//'.h5', &
                                                                      this % myRank, this % nProcs ) 

         PRINT*, '  Loading Operator : '//TRIM(fileBase)//'/diffusion.'//fileIDChar//'.h5'
         CALL this % solution % diffusionOp % ReadCRSMatrix_HDF5( TRIM(fileBase)//'/diffusion'//fileIDChar//'.h5', &
                                                                      this % myRank, this % nProcs ) 

!         IF( this % params % waterMassTagging )THEN
!
!            WRITE(fileIDChar, '(I5.5)' ) operatorPeriod
!            oceanStateFile = TRIM(this % params % regionalOperatorDirectory)//'Ocean.'//fileIDChar//'.nc'
!            PRINT*, 'Loading new ocean state : ', TRIM(oceanStateFile)
!            CALL this % nativeSol % LoadOceanState( this % mesh, TRIM(oceanStateFile) )
!
!            ! Set any prescribed cells here
!#ifndef HAVE_MPI
!            DO iMask = 1, this % feotsMap % nMasks
!               DO iLayer = 1, this % params % nLayers
!                  
!                  iTracer = iLayer + (iMask-1)*( this % params % nLayers )
!#else
!                  iMask   = (myRank-1)/(this % params % nLayers)+1
!                  iLayer  = myRank  - (iMask-1)*this % params % nLayers
!                  iTracer = 1
!#endif
!
!                  DO m = 1, this % feotsMap % bMap(iMask) % nPCells
!
!                     dof = this % feotsMap % bMap(iMask) % prescribedCells(m)
!                     i   = this % feotsMap % dofToLocalIJK(1,dof)
!                     j   = this % feotsMap % dofToLocalIJK(2,dof)
!                     k   = this % feotsMap % dofToLocalIJK(3,dof)
!
!                     trackingVar = this % nativeSol % temperature(i,j,k)*this % statemask(1) +&
!                                   this % nativeSol % salinity(i,j,k)*this % statemask(2) +&
!                                   this % nativeSol % density(i,j,k)*this % statemask(3)                 
!                  
!                     this % solution % tracers(dof,iTracer) = 0.0_prec ! Reset prescribed values
!
!                     IF( trackingVar >= this % stateLowerBound(iLayer) .AND. &
!                         trackingVar < this % stateUpperBound(iLayer) )THEN
!                        this % solution % tracers(dof,iTracer) = 1.0_prec
!                     ENDIF
!
!                  ENDDO
!#ifndef HAVE_MPI
!               ENDDO 
!            ENDDO
!#endif
!         ENDIF
      ELSE   
         
         IF( this % params % Regional ) THEN
            CALL this % solution % CheckForNewOperator( tn, &
                                TRIM( this % params % regionalOperatorDirectory), &
                                operatorsSwapped, this % myRank, this % nProcs )
         ELSE
            CALL this % solution % CheckForNewOperator( tn, &
                                TRIM( this % params % dbRoot )//'/ops/', &
                                operatorsSwapped, this % myRank, this % nProcs )
         ENDIF

!         IF( this % params % waterMassTagging .AND. operatorsSwapped )THEN
!
!            WRITE(fileIDChar, '(I5.5)' ) this % solution % currentPeriod + 1 
!            oceanStateFile = TRIM(this % params % regionalOperatorDirectory)//'Ocean.'//fileIDChar//'.nc'
!            PRINT*, 'Loading new ocean state : ', TRIM(oceanStateFile)
!            CALL this % nativeSol % LoadOceanState( this % mesh, TRIM(oceanStateFile) )
!
!            ! Set any prescribed cells here
!#ifndef HAVE_MPI
!            DO iMask = 1, this % feotsMap % nMasks
!               DO iLayer = 1, this % params % nLayers
!                  
!                  iTracer = iLayer + (iMask-1)*( this % params % nLayers )
!#else
!                  iMask   = (myRank-1)/(this % params % nLayers)+1
!                  iLayer  = myRank  - (iMask-1)*this % params % nLayers
!                  iTracer = 1
!#endif
!
!                  DO m = 1, this % feotsMap % bMap(iMask) % nPCells
!
!                     dof = this % feotsMap % bMap(iMask) % prescribedCells(m)
!                     i   = this % feotsMap % dofToLocalIJK(1,dof)
!                     j   = this % feotsMap % dofToLocalIJK(2,dof)
!                     k   = this % feotsMap % dofToLocalIJK(3,dof)
!
!                     trackingVar = this % nativeSol % temperature(i,j,k)*this % statemask(1) +&
!                                   this % nativeSol % salinity(i,j,k)*this % statemask(2) +&
!                                   this % nativeSol % density(i,j,k)*this % statemask(3)                 
!                  
!                     this % solution % tracers(dof,iTracer) = 0.0_prec ! Reset prescribed values
!
!                     IF( trackingVar >= this % stateLowerBound(iLayer) .AND. &
!                         trackingVar < this % stateUpperBound(iLayer) )THEN
!                        this % solution % tracers(dof,iTracer) = 1.0_prec
!                     ENDIF
!
!                  ENDDO
!#ifndef HAVE_MPI
!               ENDDO 
!            ENDDO
!#endif
!         ENDIF

      ENDIF

      !DO row = 1, this % feotsMap % nCells
      !  i = this % feotsMap % dofToLocalIJK(1,row)
      !  j = this % feotsMap % dofToLocalIJK(2,row)
      !  k = this % feotsMap % dofToLocalIJK(3,row)
      !  IF( j >= 105 )THEN
      !     DO i = this % solution % transportOp % rowBounds(1,row), this % solution % transportOp % rowBounds(2,row)
      !       col = this % solution % transportOp % col(i)
      !       ii = this % feotsMap % dofToLocalIJK(1,col)
      !       jj = this % feotsMap % dofToLocalIJK(2,col)
      !       kk = this % feotsMap % dofToLocalIJK(3,col)
      !       PRINT*, 'Transport Mat Val :', row, col,&
      !                                      i, j, k, ii, jj, kk, &
      !                                      this % solution % transportOp % A(i)
    
      !     ENDDO 
      !  ENDIF
      !ENDDO
      !STOP 
      !$OMP END MASTER
      !$OMP BARRIER

 END SUBROUTINE LoadNewStates_POP_FEOTS
!
 SUBROUTINE ForwardStep_POP_FEOTS( this, tn, nTimeSteps, myRank, nProcs )
 ! S/R ForwardStep
 !
 !  Identical to "CycleIntegration", except that the operators are loaded in
 !  as a function of time
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( POP_FEOTS ), INTENT(inout) :: this
   REAL(prec), INTENT(inout)         :: tn
   INTEGER, INTENT(in)               :: nTimeSteps, myRank, nProcs


      IF( this % params % timeStepScheme == Euler )THEN

         CALL this % ForwardStepEuler( tn, nTimeSteps, myRank, nProcs )

      ELSEIF( this % params % timeStepScheme == AB2 )THEN

         CALL this % ForwardStepAB2( tn, nTimeSteps, myRank, nProcs )

      ELSEIF( this % params % timeStepScheme == AB3 )THEN

         CALL this % ForwardStepAB3( tn, nTimeSteps, myRank, nProcs )

      ENDIF


 END SUBROUTINE ForwardStep_POP_FEOTS
!
 FUNCTION VerticalMixingAction( this, x ) RESULT(Ax)
   IMPLICIT NONE
   CLASS( POP_FEOTS ) :: this
   REAL(prec) :: x(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: Ax(1:this % solution % nDOF, 1:this % solution % nTracers)
   ! Local
   INTEGER :: i, m
   REAL(prec) :: Dx(1:this % solution % nDOF, 1:this % solution % nTracers)

     Dx = DiffusiveAction( this % solution % diffusionOp, x, this % solution % nTracers, this % solution % nDOF )
     DO m = 1, this % solution % nTracers
       !$OMP DO
       DO i = 1, this % solution % nDOF
         Ax(i,m) = (1.0_prec + this % solution % volume(i))*x(i,m) - this % params % dt*Dx(i,m)
       ENDDO
       !$OMP ENDDO
     ENDDO
     

 END FUNCTION VerticalMixingAction 
    
 FUNCTION GetMax( this, x ) RESULT( maxX )
   IMPLICIT NONE
   CLASS( POP_FEOTS ) :: this
   REAL(prec) :: x(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: maxX
   ! Local
   INTEGER :: i, m

     !$OMP MASTER
     maxX = 0.0_prec 
     DO m = 1, this % solution % nTracers
       DO i = 1, this % solution % nDOF
          maxX = MAX(ABS(x(i,m)), maxX)
       ENDDO
     ENDDO
     !$OMP END MASTER

 END FUNCTION GetMax

 FUNCTION Mask( this, x ) RESULT( y )
   IMPLICIT NONE
   CLASS( POP_FEOTS ) :: this
   REAL(prec) :: x(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: y(1:this % solution % nDOF, 1:this % solution % nTracers)
   ! Local
   INTEGER :: i, m

     DO m = 1, this % solution % nTracers
       !$OMP DO
       DO i = 1, this % solution % nDOF
          y(i,m) = this % solution % mask(i,m)*x(i,m) 
       ENDDO
       !$OMP ENDDO
     ENDDO

 END FUNCTION Mask

 SUBROUTINE VerticalMixing_POP_FEOTS( this, rhs )
   ! Preconditioned Conjugate Gradient used to solve vertical mixing
   ! Algorithm taken from page 3 of
   ! http://www.cse.psu.edu/~b58/cse456/lecture20.pdf
   CLASS( POP_FEOTS ), INTENT(inout) :: this
   REAL(prec), INTENT(in)            :: rhs(1:this % solution % nDOF, 1:this % solution % nTracers)
   ! Local
   REAL(prec) :: Ax(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: r(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: z(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: rk(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: zk(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: w(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: p(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: x(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: alpha, beta, residual_magnitude, sol_magnitude, update_magnitude
   INTEGER :: m, i, iter
   REAL(prec) :: mag_divisor
  

     !$OMP BARRIER
     DO m = 1, this % solution % nTracers
       !$OMP DO
       DO i = 1, this % solution % nDOF
         x(i,m) = this % solution % tracers(i,m)
       ENDDO
       !$OMP ENDDO
     ENDDO
    
     !$OMP BARRIER
     Ax = this % VerticalMixingAction( x )
     DO m = 1, this % solution % nTracers
       !$OMP DO
       DO i = 1, this % solution % nDOF
         r(i,m)  = rhs(i,m) - Ax(i,m)
       ENDDO
       !$OMP ENDDO
     ENDDO

     r = this % Mask( r )

     DO m = 1, this % solution % nTracers
       !$OMP DO
       DO i = 1, this % solution % nDOF
         ! Invert the preconditioner
         z(i,m) = r(i,m)
         p(i,m) = z(i,m)
       ENDDO
       !$OMP ENDDO
     ENDDO  

     !$OMP BARRIER
     w = this % VerticalMixingAction( p ) 
     !$OMP BARRIER

     alpha = this % DotProduct( r, z )/this % DotProduct( p, w )

     DO m = 1, this % solution % nTracers
       !$OMP DO
       DO i = 1, this % solution % nDOF
         x(i,m) = x(i,m) + alpha*p(i,m)
         rk(i,m) = r(i,m) - alpha*w(i,m)
       ENDDO
       !$OMP ENDDO
     ENDDO
     !$OMP BARRIER
     rk = this % Mask( rk )
     residual_magnitude = this % GetMax(rk)
     !$OMP FLUSH(residual_magnitude)

     IF( residual_magnitude <= cg_tolerance )THEN
       PRINT*, 'Vertical Mixing : Final Residual :', residual_magnitude
       RETURN
     ENDIF

     DO iter = 1, cg_itermax

       DO m = 1, this % solution % nTracers
         !$OMP DO
         DO i = 1, this % solution % nDOF
           ! Invert the preconditioner
           zk(i,m) = rk(i,m)
         ENDDO
         !$OMP ENDDO
       ENDDO  

       beta = this % DotProduct(rk,zk)/this % DotProduct(r,z)
       
       DO m = 1, this % solution % nTracers
         !$OMP DO
         DO i = 1, this % solution % nDOF
           p(i,m) = zk(i,m) + beta*p(i,m)
         ENDDO
         !$OMP ENDDO
       ENDDO  

       !$OMP BARRIER
       Ax = DiffusiveAction( this % solution % diffusionOp, p, this % solution % nTracers, this % solution % nDOF )
       DO m = 1, this % solution % nTracers
         !$OMP DO
         DO i = 1, this % solution % nDOF
           w(i,m) = this % solution % mask(i,m)*( (1.0_prec + this % solution % volume(i))*p(i,m) - this % params % dt*Ax(i,m) )
         ENDDO
         !$OMP ENDDO
       ENDDO

       alpha = this % DotProduct(rk,zk)/this % DotProduct(p,w)
       p = this % Mask( p )

       DO m = 1, this % solution % nTracers
         !$OMP DO
         DO i = 1, this % solution % nDOF

           r(i,m) = rk(i,m)
           z(i,m) = zk(i,m)

           this % solution % tracers(i,m) = this % solution % tracers(i,m) + alpha*p(i,m)
           rk(i,m) = r(i,m) - alpha*w(i,m)

         ENDDO
         !$OMP ENDDO
       ENDDO
     
       !$OMP BARRIER
       rk = this % Mask( rk )
       residual_magnitude = this % GetMax(rk)
       PRINT*, 'Residual :', residual_magnitude 
       !$OMP FLUSH(residual_magnitude)
       IF( residual_magnitude <= cg_tolerance )THEN
         PRINT*, 'Vertical Mixing : Final Residual :', residual_magnitude
         RETURN
       ENDIF

     ENDDO

     PRINT*, 'POP_FEOTS_Class.F90 : VerticalMixing : Failed to converge. Final residual : ', residual_magnitude
     

 END SUBROUTINE VerticalMixing_POP_FEOTS
!
 SUBROUTINE StepForward( this, dCdt, dVdt, nTimeSteps ) 
   IMPLICIT NONE
   CLASS( POP_FEOTS ), INTENT(inout) :: this
   REAL(prec), INTENT(in)            :: dCdt(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec), INTENT(in)            :: dVdt(1:this % solution % nDOF)
   INTEGER, INTENT(in)               :: nTimeSteps
   ! Local
   INTEGER    :: maskID, m, i
   REAL(prec) :: vol(1:this % solution % nDOF)
   REAL(prec) :: rhs(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: Dx(1:this % solution % nDOF, 1:this % solution % nTracers)


        !$OMP DO
        DO i = 1, this % solution % nDOF
           ! Calculate volume correction
           vol(i) = this % solution % volume(i) + this % params % dt*dVdt(i)
        ENDDO
        !$OMP ENDDO

        DO m = 1, this % solution % nTracers
           !$OMP DO
           DO i = 1, this % solution % nDOF

             rhs(i,m)  = ((1.0_prec+this % solution % volume(i))*this % solution % tracers(i,m) + this % params % dt*(dCdt(i,m)))

             ! Set the initial guess for the vertical diffusion
             this % solution % tracers(i,m) = rhs(i,m)/(1.0_prec+vol(i))

           ENDDO
           !$OMP ENDDO
        ENDDO

        ! Update volume
        !$OMP DO
        DO i = 1, this % solution % nDOF
           this % solution % volume(i) = vol(i)
        ENDDO
        !$OMP ENDDO

        ! Need to invert ( 1 + vol(i) - D )*c = rhs with Conjugate Gradient.
        !CALL this % VerticalMixing( rhs )

        !$OMP BARRIER

 END SUBROUTINE StepForward
!
 SUBROUTINE ForwardStepEuler_POP_FEOTS( this, tn, nTimeSteps, myRank, nProcs )
 ! S/R ForwardStepEuler
 !
 !  Identical to "CycleIntegrationEuler", except that the operators are loaded in
 !  as a function of time
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( POP_FEOTS ), INTENT(inout) :: this
   REAL(prec), INTENT(inout)         :: tn
   INTEGER, INTENT(in)               :: nTimeSteps, myRank, nProcs
   ! Local
   REAL(prec) :: trm1(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: trm2(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: vol(1:this % solution % nDOF)
   REAL(prec) :: weightedTracers(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: dCdt(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: diffTendency(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: dVdt(1:this % solution % nDOF)
   INTEGER    :: iT


      DO iT = 1, nTimeSteps

         CALL this % LoadNewStates( tn, myRank )
         CALL this % solution % CalculateTendency( this % solution % tracers, tn, this % params % TracerModel, dCdt, dVdt )
         CALL this % StepForward( dCdt, dVdt, nTimeSteps )

         !$OMP MASTER
         tn = tn + this % params % dt
         !$OMP END MASTER

      ENDDO
   
 END SUBROUTINE ForwardStepEuler_POP_FEOTS
!
 SUBROUTINE ForwardStepAB2_POP_FEOTS( this, tn, nTimeSteps, myRank, nProcs )
 ! S/R ForwardStepAB2
 !
 !  Identical to "CycleIntegrationAB2", except that the operators are loaded in
 !  as a function of time
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( POP_FEOTS ), INTENT(inout) :: this
   REAL(prec), INTENT(inout)         :: tn
   INTEGER, INTENT(in)               :: nTimeSteps, myRank, nProcs
   ! Local
   REAL(prec) :: trm1(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: trm2(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: vol(1:this % solution % nDOF)
   REAL(prec) :: weightedTracers(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: dCdt(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: diffTendency(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: dVdt(1:this % solution % nDOF)
   INTEGER         :: i, j, k, m, iT, dof, iTracer, iLayer, iMask
   LOGICAL         :: operatorsSwapped
   CHARACTER(400)  :: oceanStateFile
   CHARACTER(5)    :: fileIDChar
   REAL(prec)      :: trackingVar

      DO iT = 1, nTimeSteps

         CALL this % LoadNewStates( tn, myRank )

         IF( iT == 1 )THEN! First order Euler
               !$OMP DO COLLAPSE(2)
               DO m = 1, this % solution % nTracers
                  DO i = 1, this % solution % nDOF
                     weightedTracers(i,m) = this % solution % tracers(i,m)
                     trm1(i,m) = this % solution % tracers(i,m)
                  ENDDO
               ENDDO
               !$OMP ENDDO

         ELSE ! Second Order Adams Bashforth
               !$OMP DO COLLAPSE(2)
               DO m = 1, this % solution % nTracers
                  DO i = 1, this % solution % nDOF
                     weightedTracers(i,m) = (3.0_prec*this % solution % tracers(i,m) - trm1(i,m))*0.5_prec
                     trm2(i,m) = trm1(i,m)
                     trm1(i,m) = this % solution % tracers(i,m)
                  ENDDO
               ENDDO
               !$OMP ENDDO
         ENDIF

         !$OMP BARRIER
         CALL this % solution % CalculateTendency( weightedTracers, tn, this % params % TracerModel, dCdt, dVdt)
         !$OMP BARRIER

         CALL this % StepForward( dCdt, dVdt, nTimeSteps )

         !$OMP MASTER
         tn = tn + this % params % dt
         !$OMP END MASTER

      ENDDO
   
 END SUBROUTINE ForwardStepAB2_POP_FEOTS
!
 SUBROUTINE ForwardStepAB3_POP_FEOTS( this, tn, nTimeSteps, myRank, nProcs )
 ! S/R ForwardStepAB3
 !
 !  Identical to "CycleIntegrationAB3", except that the operators are loaded in
 !  as a function of time
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( POP_FEOTS ), INTENT(inout) :: this
   REAL(prec), INTENT(inout)         :: tn
   INTEGER, INTENT(in)               :: nTimeSteps, myRank, nProcs
   ! Local
   REAL(prec) :: trm1(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: trm2(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: vol(1:this % solution % nDOF)
   REAL(prec) :: weightedTracers(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: dCdt(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: diffTendency(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: dVdt(1:this % solution % nDOF)
   INTEGER         :: i, j, k, m, iT, dof, iTracer, iLayer, iMask
   LOGICAL         :: operatorsSwapped
   CHARACTER(400)  :: oceanStateFile
   CHARACTER(5)    :: fileIDChar
   REAL(prec)      :: trackingVar

      DO iT = 1, nTimeSteps

         CALL this % LoadNewStates( tn, myRank )


         IF( iT == 1 )THEN! First order Euler
               !$OMP DO COLLAPSE(2)
               DO m = 1, this % solution % nTracers
                  DO i = 1, this % solution % nDOF
                     weightedTracers(i,m) = this % solution % tracers(i,m)
                     trm1(i,m) = this % solution % tracers(i,m)
                  ENDDO
               ENDDO
               !$OMP ENDDO

         ELSEIF( iT == 2 )THEN ! Second Order Adams Bashforth
               !$OMP DO COLLAPSE(2)
               DO m = 1, this % solution % nTracers
                  DO i = 1, this % solution % nDOF
                     weightedTracers(i,m) = (3.0_prec*this % solution % tracers(i,m) - trm1(i,m))*0.5_prec
                     trm2(i,m) = trm1(i,m)
                     trm1(i,m) = this % solution % tracers(i,m)
                  ENDDO
               ENDDO
               !$OMP ENDDO

         ELSE! Third Order Adams Bashforth
               !$OMP DO COLLAPSE(2)
               DO m = 1, this % solution % nTracers
                  DO i = 1, this % solution % nDOF
                     weightedTracers(i,m) = (23.0_prec*this % solution % tracers(i,m) - 16.0_prec*trm1(i,m) + 5.0_prec*trm2(i,m))/12.0_prec
                     trm2(i,m) = trm1(i,m)
                     trm1(i,m) = this % solution % tracers(i,m)
                  ENDDO
               ENDDO
               !$OMP ENDDO
         ENDIF

         !$OMP BARRIER
         CALL this % solution % CalculateTendency( weightedTracers, tn, this % params % TracerModel, dCdt, dVdt)
         !$OMP BARRIER

         CALL this % StepForward( dCdt, dVdt, nTimeSteps )

         !$OMP MASTER
         tn = tn + this % params % dt
         !$OMP END MASTER

      ENDDO
   
 END SUBROUTINE ForwardStepAB3_POP_FEOTS
!
 SUBROUTINE CycleIntegration_POP_FEOTS( this, myRank )
 ! S/R ForwardStep
 !
 !  Identical to "CycleIntegration", except that the operators are loaded in
 !  as a function of time
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( POP_FEOTS ), INTENT(inout) :: this
   INTEGER, INTENT(in)               :: myRank

      IF( this % params % timeStepScheme == Euler )THEN

         CALL this % CycleIntegrationEuler( myRank )

      ELSEIF( this % params % timeStepScheme == AB2 )THEN

         CALL this % CycleIntegrationAB2( myRank )

      ELSEIF( this % params % timeStepScheme == AB3 )THEN

         CALL this % CycleIntegrationAB3( myRank )

      ENDIF

 END SUBROUTINE CycleIntegration_POP_FEOTS
!
 SUBROUTINE CycleIntegrationEuler_POP_FEOTS( this, myRank )
 ! S/R CycleIntegrationEuler
 !
 !    This subroutine integrates the tracer system from t=0 to t=nPeriods*opPeriod, where 
 !    "nPeriods" is the number of transport operators we have (in time) and "opPeriod" is the period
 !    of time associated with each operator.
 !    For example, if we have 12 operators that each have a 1-month period, the system is integrated
 !    over 1 year.
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( POP_FEOTS ), INTENT(inout) :: this
   INTEGER, INTENT(in)               :: myRank
   ! Local
   REAL(prec) :: trm1(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: trm2(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: vol(1:this % solution % nDOF)
   REAL(prec) :: weightedTracers(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: dCdt(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: diffTendency(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: dVdt(1:this % solution % nDOF)
   REAL(prec) :: tn, dt
   INTEGER    :: nPeriods, nSteps, i, j, k, m 
   CHARACTER(5) :: periodChar
   CHARACTER(500) :: fileBase

      nSteps = INT( this % solution % opPeriod/ this % solution % dt )

      DO j = 1, this % params % nOperatorsPerCycle

         CALL this % LoadNewStates( tn, myRank, j )

         DO i = 1, nSteps

            tn = REAL(i-1,prec)*dt

            !$OMP BARRIER
            CALL this % solution % CalculateTendency( this % solution % tracers, tn, this % params % TracerModel, dCdt, dVdt)
            !$OMP BARRIER
         
            CALL this % StepForward( dCdt, dVdt, 1 )

         ENDDO

         ! This last step ensures we end on t=tCycle
         tn = REAL(nSteps,prec)*this % params % dt
         dt = this % params % dt
         this % params % dt = this % params % operatorPeriod -tn

         CALL this % solution % CalculateTendency( this % solution % tracers, tn, this % params % TracerModel, dCdt, dVdt )

         CALL this % StepForward( dCdt, dVdt, 1 )

         this % params % dt = dt

      ENDDO
 
 END SUBROUTINE CycleIntegrationEuler_POP_FEOTS
!
 SUBROUTINE CycleIntegrationAB2_POP_FEOTS( this, myRank )
 ! S/R CycleIntegrationAB2
 !
 !    This subroutine integrates the tracer system from t=0 to t=nPeriods*opPeriod, where 
 !    "nPeriods" is the number of transport operators we have (in time) and "opPeriod" is the period
 !    of time associated with each operator.
 !    For example, if we have 12 operators that each have a 1-month period, the system is integrated
 !    over 1 year.
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( POP_FEOTS ), INTENT(inout) :: this
   INTEGER, INTENT(in)               :: myRank
   ! Local
   REAL(prec) :: trm1(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: trm2(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: vol(1:this % solution % nDOF)
   REAL(prec) :: weightedTracers(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: dCdt(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: diffTendency(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: dVdt(1:this % solution % nDOF)
   REAL(prec) :: tn, dt
   INTEGER    :: nPeriods, nSteps, i, j, k, m 
   CHARACTER(5) :: periodChar
   CHARACTER(500) :: fileBase

      nSteps = INT( this % solution % opPeriod/ this % solution % dt )

      DO j = 1, this % params % nOperatorsPerCycle

         CALL this % LoadNewStates( tn, myRank, j )

         DO i = 1, nSteps

            tn = REAL(i-1,prec)*dt

            IF( i == 1 )THEN
               !$OMP DO COLLAPSE(2)
               DO m = 1, this % solution % nTracers
                  DO k = 1, this % solution % nDOF
                     weightedTracers(k,m) = this % solution % tracers(k,m)
                     trm1(k,m) = this % solution % tracers(k,m)
                  ENDDO
               ENDDO
               !$OMP ENDDO
            ELSE
               !$OMP DO COLLAPSE(2)
               DO m = 1, this % solution % nTracers
                  DO k = 1, this % solution % nDOF
                     weightedTracers(k,m) = (3.0_prec*this % solution % tracers(k,m) - trm1(k,m))*0.5_prec
                     trm2(k,m) = trm1(k,m)
                     trm1(k,m) = this % solution % tracers(k,m)
                  ENDDO
               ENDDO
               !$OMP ENDDO
            ENDIF

            !$OMP BARRIER
            CALL this % solution % CalculateTendency( weightedTracers, tn, this % params % TracerModel, dCdt, dVdt)
            !$OMP BARRIER

            CALL this % StepForward( dCdt, dVdt, 1 )
            
         ENDDO

         ! This last step ensures we end on t=tCycle
         tn = REAL(nSteps,prec)*this % params % dt
         dt = this % params % dt
         this % params % dt = this % params % operatorPeriod -tn

         CALL this % solution % CalculateTendency( this % solution % tracers, tn, this % params % TracerModel, dCdt, dVdt)

         CALL this % StepForward( dCdt, dVdt, 1 )

         this % params % dt = dt
      ENDDO
 
 END SUBROUTINE CycleIntegrationAB2_POP_FEOTS
!
 SUBROUTINE CycleIntegrationAB3_POP_FEOTS( this, myRank )
 ! S/R CycleIntegrationAB3
 !
 !    This subroutine integrates the tracer system from t=0 to t=nPeriods*opPeriod, where 
 !    "nPeriods" is the number of transport operators we have (in time) and "opPeriod" is the period
 !    of time associated with each operator.
 !    For example, if we have 12 operators that each have a 1-month period, the system is integrated
 !    over 1 year.
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( POP_FEOTS ), INTENT(inout) :: this
   INTEGER, INTENT(in)               :: myRank
   ! Local
   REAL(prec) :: trm1(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: trm2(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: vol(1:this % solution % nDOF)
   REAL(prec) :: weightedTracers(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: dCdt(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: diffTendency(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: dVdt(1:this % solution % nDOF)
   REAL(prec) :: tn, dt
   INTEGER    :: nPeriods, nSteps, i, j, k, m 
   CHARACTER(5) :: periodChar
   CHARACTER(500) :: fileBase

      nSteps = INT( this % solution % opPeriod/ this % solution % dt )

      DO j = 1, this % params % nOperatorsPerCycle

         CALL this % LoadNewStates( tn, myRank, j )

         DO i = 1, nSteps

            tn = REAL(i-1,prec)*dt

            IF( i == 1 )THEN ! First order Euler
               !$OMP DO COLLAPSE(2)
               DO m = 1, this % solution % nTracers
                  DO k = 1, this % solution % nDOF
                     weightedTracers(k,m) = this % solution % tracers(k,m)
                     trm1(k,m) = this % solution % tracers(k,m)
                  ENDDO
               ENDDO
               !$OMP ENDDO
            ELSEIF( i == 2 )THEN
               !$OMP DO COLLAPSE(2)
               DO m = 1, this % solution % nTracers
                  DO k = 1, this % solution % nDOF
                     weightedTracers(k,m) = (3.0_prec*this % solution % tracers(k,m) - trm1(k,m))*0.5_prec
                     trm2(k,m) = trm1(k,m)
                     trm1(k,m) = this % solution % tracers(k,m)
                  ENDDO
               ENDDO
               !$OMP ENDDO
            ELSE
               !$OMP DO COLLAPSE(2)
               DO m = 1, this % solution % nTracers
                  DO k = 1, this % solution % nDOF
                     weightedTracers(k,m) = (23.0_prec*this % solution % tracers(k,m) - 16.0_prec*trm1(k,m) + 5.0_prec*trm2(k,m))/12.0_prec
                     trm2(k,m) = trm1(k,m)
                     trm1(k,m) = this % solution % tracers(k,m)
                  ENDDO
               ENDDO
               !$OMP ENDDO
            ENDIF

            !$OMP BARRIER
            CALL this % solution % CalculateTendency( weightedTracers, tn, this % params % TracerModel, dCdt, dVdt )
            !$OMP BARRIER
            
            CALL this % StepForward( dCdt, dVdt, 1 )
            
         ENDDO

         ! This last step ensures we end on t=tCycle
         tn = REAL(nSteps,prec)*this % params % dt
         dt = this % params % dt
         this % params % dt = this % params % operatorPeriod -tn

         CALL this % solution % CalculateTendency( this % solution % tracers, tn, this % params % TracerModel, dCdt, dVdt)

         CALL this % StepForward( dCdt, dVdt, 1 )

         this % params % dt = dt

      ENDDO
 
 END SUBROUTINE CycleIntegrationAB3_POP_FEOTS
!
 FUNCTION DotProduct_POP_FEOTS( this, x, y ) RESULT( xdoty )
 !
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( POP_FEOTS ) :: this
   REAL(prec)        :: x(1:this % solution % nDOF,1:this % solution % nTracers)
   REAL(prec)        :: y(1:this % solution % nDOF,1:this % solution % nTracers)
   REAL(prec)        :: xdoty
   ! Local
   INTEGER    :: i, j

      !$OMP MASTER
      xdoty = 0.0_prec
      DO j = 1, this % solution % nTracers
         DO i = 1, this % solution % nDOF
            xdoty = xdoty + x(i,j)*y(i,j)*this % solution % mask(i,j)
         ENDDO
      ENDDO
      !$OMP END MASTER

 END FUNCTION DotProduct_POP_FEOTS
!
! ================================================================================================ !
!                                  Equilibration Routines
! ================================================================================================ !
!
 FUNCTION JacobianAction_POP_FEOTS( this, x, Gx, v, myRank ) RESULT( Jv )
 !  JacobianAction
 !
 !    This subroutine calculates the action of the Jacobian matrix for the fixed point map
 !       G(x) = M(x) - x, 
 !    where M(x) takes in a solution "x" an integrates the system through all periods once - 
 !    - this is a call to "CycleIntegration". The Jacobian action is approximated using
 !    J(v) = (G(x + e*v) - G(x))/e, where "e" is a given constant in the parameters file
 !    e ~~> params % JacobianStepSize
 !
 !    On input, this % solution corresponds to x
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( POP_FEOTS ) :: this
   REAL(prec)        :: x(1:this % solution % nDOF,1:this % solution % nTracers)
   REAL(prec)        :: Gx(1:this % solution % nDOF,1:this % solution % nTracers)
   REAL(prec)        :: v(1:this % solution % nDOF,1:this % solution % nTracers)
   REAL(prec)        :: Jv(1:this % solution % nDOF,1:this % solution % nTracers)
   INTEGER           :: myRank
   ! Local
   INTEGER           :: i, j
   REAL(prec)        :: scFac, vmag, e

      scFac = 0.0_prec
      vmag  = 0.0_prec
      !$OMP DO COLLAPSE(2) REDUCTION(+:scFac,vmag)
      DO j = 1, this % solution % nTracers
         DO i = 1, this % solution % nDOF 
            scFac = scFac + this % params % JacobianStepSize*( abs(x(i,j)) )
            vmag  = vmag + v(i,j)**2
         ENDDO
      ENDDO
      !$OMP ENDDO

      !$OMP BARRIER     

      !$OMP MASTER  
      scFac = scFac + this % params % JacobianStepSize
      e = scFac/( REAL( this % solution % nDOF*this % solution % nTracers, prec)*SQRT( vmag ) )
      !$OMP END MASTER

      !$OMP BARRIER     

      !$OMP DO COLLAPSE(2)      
      DO j = 1, this % solution % nTracers
         DO i = 1, this % solution % nDOF 
            this % solution % tracers(i,j) = x(i,j) + e*v(i,j)
         ENDDO
      ENDDO
      !$OMP ENDDO
      !$OMP FLUSH( this )

      CALL this % CycleIntegrationAB3( myRank )
      !$OMP BARRIER

      !$OMP DO COLLAPSE(2)      
      DO j = 1, this % solution % nTracers
         DO i = 1, this % solution % nDOF 
            Jv(i,j) = ( this % solution % tracers(i,j) - (x(i,j) + e*v(i,j)) - Gx(i,j) )/e
         ENDDO
      ENDDO
      !$OMP ENDDO
      !$OMP FLUSH( Jv )

 END FUNCTION JacobianAction_POP_FEOTS
!
 SUBROUTINE SolveGMRES_POP_FEOTS( this, x, Gx, dx, resi, ioerr, myRank )
 !  S/R SolveGMRES
 !
 !  This subroutine solves the system J(dx) = Gx using the preconditioned GMRES.
 !  The matrix action, residual, and preconditioning routines are supplied by a non-abstracted 
 !  type-extension of PCGMRES. These routines should return an array indexed from 1 to 
 !  nDOF. Thus, in addition to a MatrixAction and Residual, the user should map their data-structure
 !  to a 1-D array so that communication to this routine is possible.  
 !
 !  The flavor of GMRES implemented here is the "Flexible GMRES with variable right preconditioning".
 !  See pp.74-76 of "Iterative Krylov Methods for Large Linear Systems" by Henk A. van der Vorst
 !
 !  On output ioerr is set to an error checking flag. 
 !  If ioerr ==  0, the method converged within the maximum number of iterations.
 !     ioerr == -1, the method did not converge within the maximum number of iterations.
 !     ioerr == -2, something that is not caught by the current construct happened.
 !  
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( POP_FEOTS ), INTENT(inout) :: this
   REAL(prec), INTENT(in)           :: x(1:this % solution % nDOF,1:this % solution % nTracers)
   REAL(prec), INTENT(in)           :: Gx(1:this % solution % nDOF,1:this % solution % nTracers)
   REAL(prec), INTENT(out)          :: dx(1:this % solution % nDOF,1:this % solution % nTracers)
   REAL(prec), INTENT(out)          :: resi(0:this % params % nResi)
   INTEGER, INTENT(out)             :: ioerr
   INTEGER, INTENT(in)              :: myRank
   ! LOCAL
   INTEGER    :: i, j, k, l
   INTEGER    :: ii, jj
   INTEGER    :: nIt, m, nr 
   REAL(prec) :: TOL
   REAL(prec) :: r(1:this % solution % nDOF, &
                   1:this % solution % nTracers)
   REAL(prec) :: v(1:this % solution % nDOF, &
                   1:this % solution % nTracers, &
                   1:this % params % mInnerItersGMRES+1)
   REAL(prec) :: w(1:this % solution % nDOF, &
                   1:this % solution % nTracers)
   REAL(prec) :: z(1:this % solution % nDOF, &
                   1:this % solution % nTracers, &
                   1:this % params % mInnerItersGMRES+1)
   REAL(prec) :: localz(1:this % solution % nDOF, &
                        1:this % params % mInnerItersGMRES)
   REAL(prec) :: rho(1:this % params % mInnerItersGMRES, &
                     1:this % params % mInnerItersGMRES)
   REAL(prec) :: h(1:this % params % mInnerItersGMRES+1, &
                   1:this % params % mInnerItersGMRES)
   REAL(prec) :: c(1:this % params % mInnerItersGMRES)
   REAL(prec) :: s(1:this % params % mInnerItersGMRES)
   REAL(prec) :: bhat(1:this % params % mInnerItersGMRES+1)
   REAL(prec) :: y(1:this % params % mInnerItersGMRES+1)
   REAL(prec) :: b, d, g, r0, rc, hki
   REAL(prec) :: e, vmag, scFac

      ioerr = -2
      nIt = this % params % maxItersGMRES
      m   = this % params % mInnerItersGMRES
      TOL = this % params % toleranceGMRES
      
      ! Assume that the initial guess is dx = 0
      r = -Gx
      r0 = sqrt( this % DotProduct( r, r ) )
     
      l = 0 
      resi = 0.0_prec
      resi(l) = r0
      
      dx   = 0.0_prec
      v    = 0.0_prec
      rho  = 0.0_prec
      bhat = 0.0_prec 
      s    = 0.0_prec
      c    = 0.0_prec
      y    = 0.0_prec
   
      !$OMP PARALLEL
      DO j = 1,nIt
          b = 0.0_prec
         !$OMP DO COLLAPSE(2) REDUCTION(+:b)
         DO jj = 1, this % solution % nTracers
            DO ii = 1, this % solution % nDOF
               b = b + r(ii,jj)**2
            ENDDO
         ENDDO
         !$OMP ENDDO 
 
         !$OMP BARRIER
         !$OMP MASTER
         b = sqrt( b )
         !$OMP END MASTER
          
         !$OMP DO COLLAPSE(2)
         DO jj = 1, this % solution % nTracers
            DO ii = 1, this % solution % nDOF
               v(ii,jj,1) = r(ii,jj)/b
            ENDDO
         ENDDO
         !$OMP ENDDO
         !$OMP FLUSH( v )

         bhat(1)  = b

         DO i = 1, m
            !$OMP MASTER
            l = l+1
            nr = i
            !$OMP END MASTER

            ! Applying the precondtioner (right now it is the identity, ie, no preconditioning)
            !$OMP DO COLLAPSE(2)
            DO jj = 1, this % solution % nTracers
               DO ii = 1, this % solution % nDOF
                  z(ii,jj,i) = v(ii,jj,i)
               ENDDO
            ENDDO
            !$OMP ENDDO 
            !$OMP FLUSH( z )

            !z(:,i) = this % PreconditionInvert( v(:,i) )

            ! The first step in GMRES is to build the orthogonal basis up to order "i"
            ! with the accompanying upper hessenburg matrix.
!            w = this % JacobianAction( x, Gx, z(:,:,i) )
            !****************************** Manually Inlined "JacobianAction" ***************************** !
            scFac = 0.0_prec
            vmag  = 0.0_prec
            !$OMP DO COLLAPSE(2) REDUCTION(+:scFac,vmag)
            DO jj = 1, this % solution % nTracers
               DO ii = 1, this % solution % nDOF 
                  scFac = scFac + this % params % JacobianStepSize*( abs(x(ii,jj)) )
                  vmag  = vmag + z(ii,jj,i)**2
               ENDDO
            ENDDO
            !$OMP ENDDO

            !$OMP BARRIER     

            !$OMP MASTER  
            scFac = scFac + this % params % JacobianStepSize
            e = scFac/( REAL( this % solution % nDOF*this % solution % nTracers, prec)*SQRT( vmag ) )
            !$OMP END MASTER

            !$OMP BARRIER     

            !$OMP DO COLLAPSE(2)      
            DO jj = 1, this % solution % nTracers
               DO ii = 1, this % solution % nDOF 
                  this % solution % tracers(ii,jj) = x(ii,jj) + e*z(ii,jj,i)
               ENDDO
            ENDDO
            !$OMP ENDDO
            !$OMP FLUSH( this )

            CALL this % CycleIntegration( myRank )
            !$OMP BARRIER

            !$OMP DO COLLAPSE(2)      
            DO jj = 1, this % solution % nTracers
               DO ii = 1, this % solution % nDOF 
                  w(ii,jj) = ( this % solution % tracers(ii,jj) - (x(ii,jj) + e*z(ii,jj,i)) - Gx(ii,jj) )/e
               ENDDO
            ENDDO
            !$OMP ENDDO
            !$OMP FLUSH( w )
            !******************************            END                  " ***************************** !
            !****************************** Manually Inlined "JacobianAction" ***************************** !


            ! The new basis vector is obtained by multiplying the previous basis vector by the matrix
            ! and orthogonalizing wrt to all of the previous basis vectors using a Gram-Schmidt process.
            DO k = 1, i
                              
               hki = 0.0_prec

               !$OMP DO COLLAPSE(2) REDUCTION(+:hki)
               DO jj = 1, this % solution % nTracers
                  DO ii = 1, this % solution % nDOF
                     hki = hki + v(ii,jj,k)*w(ii,jj)
                  ENDDO
               ENDDO
               !$OMP ENDDO

               h(k,i) = hki 
               !$OMP DO COLLAPSE(2)
               DO jj = 1, this % solution % nTracers
                  DO ii = 1, this % solution % nDOF
                     w(ii,jj) = w(ii,jj) - h(k,i)*v(ii,jj,k)
                  ENDDO
               ENDDO
               !$OMP ENDDO
               !$OMP FLUSH( w )

            ENDDO

                              
            hki = 0.0_prec
            !$OMP DO COLLAPSE(2) REDUCTION(+:hki)
            DO jj = 1, this % solution % nTracers
               DO ii = 1, this % solution % nDOF
                  hki = hki + w(ii,jj)**2
               ENDDO
            ENDDO
            !$OMP ENDDO
            h(i+1,i) = sqrt(hki)

            IF( AlmostEqual( h(i+1,i), 0.0_prec )  )THEN
               EXIT
            ENDIF

            !$OMP DO COLLAPSE(2)
            DO jj = 1, this % solution % nTracers
               DO ii = 1, this % solution % nDOF
                  v(ii,jj,i+1) = w(ii,jj)/h(i+1,i)
               ENDDO
            ENDDO
            !$OMP ENDDO
            !$OMP FLUSH(v)

            !$OMP MASTER
            rho(1,i) = h(1,i)
     
            ! Givens rotations are applied to the upper hessenburg matrix and to the residual vectors
            ! that are formed from the orthogonalization process. Here, they are done "on-the-fly"
            ! as opposed to building the entire upper hessenburg matrix and orthonormal basis
            ! before performing the rotations. This way, we can also tell if we have found an exact
            ! solution ( if h(i+1,i) = 0 ) with a smaller subspace than size m.
            DO k = 2, i
               g          = c(k-1)*rho(k-1,i) + s(k-1)*h(k,i)
               rho(k,i)   = -s(k-1)*rho(k-1,i) + c(k-1)*h(k,i)
               rho(k-1,i) = g 
            ENDDO

            ! Here the coefficients of the Givens rotation matrix are computed
            d = sqrt( rho(i,i)**2 + h(i+1,i)**2 )
            c(i) = rho(i,i)/d
            s(i) = h(i+1,i)/d

            rho(i,i) = c(i)*rho(i,i) + s(i)*h(i+1,i)
            ! And applied to the residual vector
            bhat(i+1) = -s(i)*bhat(i)
            bhat(i)   = c(i)*bhat(i)

            rc = abs( bhat(i+1) )
            resi(l) = rc
            !$OMP END MASTER
            !$OMP BARRIER

            IF( rc/r0 <= TOL )THEN
               EXIT
            ENDIF
         ENDDO

         !$OMP MASTER
         IF( rc/r0 > TOL )THEN
            nr = m
         ENDIF

         ! Back-substitution of the tridiagonal matrix that resulted from the rotations
         y(nr) = bhat(nr)/rho(nr,nr)
         DO k = nr-1, 1, -1
            y(k) = bhat(k)
            DO i = k+1, nr
               y(k) = y(k) - rho(k,i)*y(i)
            ENDDO
            y(k) = y(k)/rho(k,k)
         ENDDO
         !$OMP END MASTER
  
         !$OMP BARRIER
     
         DO jj = 1, nr
            !$OMP DO COLLAPSE(2)
            DO i = 1, this % solution % nTracers
               DO ii = 1, this % solution % nDOF
                  dx(ii,i) = dx(ii,i) + z(ii,i,jj)*y(jj)
               ENDDO
            ENDDO
            !$OMP ENDDO
         ENDDO

         !$OMP BARRIER
         IF( rc/r0 <= TOL )THEN
            ioerr = l
            EXIT
         ENDIF
         !r = -this % JacobianAction( x, Gx, dx )
         !****************************** Manually Inlined "JacobianAction" ***************************** !
            scFac = 0.0_prec
            vmag  = 0.0_prec
            !$OMP DO COLLAPSE(2) REDUCTION(+:scFac,vmag)
            DO jj = 1, this % solution % nTracers
               DO ii = 1, this % solution % nDOF 
                  scFac = scFac + this % params % JacobianStepSize*( abs(x(ii,jj)) )
                  vmag  = vmag + dx(ii,jj)**2
               ENDDO
            ENDDO
            !$OMP ENDDO

            !$OMP BARRIER     

            !$OMP MASTER  
            scFac = scFac + this % params % JacobianStepSize
            e = scFac/( REAL( this % solution % nDOF*this % solution % nTracers, prec)*SQRT( vmag ) )
            !$OMP END MASTER

            !$OMP BARRIER     

            !$OMP DO COLLAPSE(2)      
            DO jj = 1, this % solution % nTracers
               DO ii = 1, this % solution % nDOF 
                  this % solution % tracers(ii,jj) = x(ii,jj) + e*dx(ii,jj)
               ENDDO
            ENDDO
            !$OMP ENDDO
            !$OMP FLUSH( this )

            CALL this % CycleIntegration( myRank )
            !$OMP BARRIER

            !$OMP DO COLLAPSE(2)      
            DO jj = 1, this % solution % nTracers
               DO ii = 1, this % solution % nDOF 
                  r(ii,jj) = ( this % solution % tracers(ii,jj) - (x(ii,jj) + e*dx(ii,jj)) - Gx(ii,jj) )/e
               ENDDO
            ENDDO
            !$OMP ENDDO
            !$OMP FLUSH( r )
         !******************************            END                  " ***************************** !
         !****************************** Manually Inlined "JacobianAction" ***************************** !

         !$OMP DO COLLAPSE(2)
         DO jj = 1, this % solution % nTracers
            DO ii = 1, this % solution % nDOF
               r(ii,jj) = r(ii,jj)-Gx(ii,jj)
            ENDDO
         ENDDO
         !$OMP ENDDO
         !$OMP FLUSH( r )
         
      ENDDO 
      !$OMP END PARALLEL

      IF( rc/r0 > TOL )THEN
         PRINT*, 'MODULE IterativeSolvers : GMRES failed to converge '
         PRINT*, 'Last L-2 residual : ', sqrt( this % DotProduct(r,r))
         ioerr=-l
      ENDIF   
    

 END SUBROUTINE SolveGMRES_POP_FEOTS
!
 SUBROUTINE JFNK_POP_FEOTS( this, myRank )
 !
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( POP_FEOTS ), INTENT(inout) :: this
   INTEGER, INTENT(in)               :: myRank
   ! Local 
   REAL(prec)   :: x(1:this % solution % nDOF,1:this % solution % nTracers)
   REAL(prec)   :: dx(1:this % solution % nDOF,1:this % solution % nTracers)
   REAL(prec)   :: Gx(1:this % solution % nDOF,1:this % solution % nTracers)
   REAL(prec)   :: residual(0:this % params % nResi)
   REAL(prec)   :: tolerance, dxMag, xMag, Gx0mag, Gxmag
   INTEGER      :: iter, iterMax, ioerr, nIt, m, i, j, k, tk, fUnit, nr
   INTEGER      :: maskID 
   CHARACTER(4) :: iterChar


      iterMax   = this % params % maxItersJFNK
      tolerance = this % params % toleranceJFNK
      nIt       = this % params % maxItersGMRES
      m         = this % params % mInnerItersGMRES
      ! Assign the initial guess
      x  = this % solution % tracers
      dx = 0.0_prec

      tk = 1

      OPEN( UNIT   = NewUnit(fUnit), &
            FILE   = 'Residual.curve', &
            FORM   = 'FORMATTED', &
            STATUS = 'REPLACE', &
            ACTION = 'WRITE' )
      WRITE( fUnit, * )'#Residual'

      CALL this % CycleIntegration( myRank )
      Gx = this % solution % tracers - x
      Gx0Mag = SQRT( this % DotProduct( Gx, Gx ) )

      DO iter = 1, iterMax
         PRINT*, 'S/R JFNK : Iterate = ', iter

         ! The first call to "CycleIntegration" calculates G(x_i),
         ! which is the RHS of the Linearized fixed point problem
         !    J*dx = -Gx

         
         PRINT*, 'S/R JFNK : Call to SolveGMRES '
         CALL this % SolveGMRES( x, Gx, dx, residual, ioerr, myRank )
         ! Check for convergence
         dxMag = SQRT( this % DotProduct( dx, dx ) )
         xMag = SQRT( this % DotProduct( x, x ) )
         ! Update the solution
         x = x + dx

         IF( ioerr < 0)THEN
            PRINT*, 'S/R JFNK : Call to SolveGMRES : Unsuccessful, ioerr = ',ioerr
         ELSE
            PRINT*, 'S/R JFNK : Call to SolveGMRES : Successful! '
         ENDIF
         nr = ABS( ioerr )
         ! Write the Residual to file
         k = 1
         
         DO i = 1, nr
            WRITE( fUnit, * ) tk, residual(i)
            tk = tk+1
         ENDDO

         this % solution % tracers = x
         CALL this % CycleIntegration( myRank )
         Gx = this % solution % tracers - x
         GxMag = SQRT( this % DotProduct( Gx, Gx ) )
         PRINT*, xMag, dxMag, GxMag

      !   IF( dxMag/xMag <= tolerance .AND. GxMag <= tolerance )THEN
         IF( GxMag <= tolerance )THEN
            PRINT*, ' Module POP_FEOTS_Class.f90 : S/R JFNK_POP_FEOTS : Solution found in ', iter, ' iterates.'
            PRINT*, ' Final Residual : ', GxMag, Gx0Mag
            PRINT*, ' Final Solution Update magnitude : ', dxMag, xMag
            EXIT
         ENDIF

         IF( MOD( iter, this % params % nStepsPerDump ) == 0 )THEN
         
            CALL this % MapTracerFromDOF( ) ! Map from DOF to native storage
            CALL this % nativeSol % InitializeForNetCDFWrite( this % params % TracerModel, &
                                                              this % mesh, &
                                                              'Tracer.pickup.nc', &
                                                              .false. )
            CALL this % nativeSol % WriteNetCDFRecord( this % mesh, 1 )
            CALL this % nativeSol % FinalizeNetCDF( )
         ENDIF
      ENDDO

      IF( GxMag <= tolerance )THEN
         ! Converged, write the solution !
         PRINT*, ' Module POP_FEOTS_Class.f90 : S/R JFNK_POP_FEOTS : Solution found in ', iter, ' iterates.'
         PRINT*, ' Final Residual : ', GxMag, Gx0Mag
         PRINT*, ' Final Solution Update magnitude : ', dxMag, xMag
         this % solution % tracers = x
         CALL this % MapTracerFromDOF( ) ! Map from DOF to native storage
         CALL this % nativeSol % InitializeForNetCDFWrite( this % params % TracerModel, &
                                                           this % mesh, &
                                                           'Tracer.equilibrium.nc', &
                                                           .false. )
         CALL this % nativeSol % WriteNetCDFRecord( this % mesh, 1 )
         CALL this % nativeSol % FinalizeNetCDF( )
         CLOSE(fUnit)
      ELSE
         ! Not converged, write a pickup file !
         !WRITE( iterChar, '(I4.4)' ) iterMax + this % params % cycleStart
         this % solution % tracers = x 
         CALL this % MapTracerFromDOF( )
         !
         CALL this % nativeSol % InitializeForNetCDFWrite( this % params % TracerModel, &
                                                           this % mesh, &
                                                           'Tracer.pickup.nc', &
                                                           .false. )
         CALL this % nativeSol % WriteNetCDFRecord( this % mesh, 1 )
         CALL this % nativeSol % FinalizeNetCDF( )


         PRINT*, ' Module POP_FEOTS_Class.f90 : S/R JFNK_POP_FEOTS : Solution not found in ', iter, ' iterates.'
         PRINT*, ' Final Residual : ', GxMag
         CLOSE(fUnit)
      ENDIF
         

 END SUBROUTINE JFNK_POP_FEOTS
!
END MODULE POP_FEOTS_Class
