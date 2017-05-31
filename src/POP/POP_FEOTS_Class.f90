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

 IMPLICIT NONE
#ifdef HAVE_MPI
INCLUDE 'mpif.h'
#endif


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


   TYPE POP_FEOTS
     
      TYPE( TracerStorage )     :: solution
      TYPE( POP_Regional )      :: regionalMaps
      TYPE( POP_Mesh )          :: mesh
      TYPE( POP_Params )        :: params
      TYPE( POP_Native )        :: nativeSol

      ! Water Mass Tagging
      REAL(prec)                :: stateMask(1:3)
      REAL(prec), ALLOCATABLE   :: stateLowerBound(:)
      REAL(prec), ALLOCATABLE   :: stateUpperBound(:)

      CONTAINS

      PROCEDURE :: Build => Build_POP_FEOTS
      PROCEDURE :: Trash => Trash_POP_FEOTS

      PROCEDURE :: SetupSettlingOperator => SetupSettlingOperator_POP_FEOTS

#ifdef HAVE_MPI
      PROCEDURE :: BroadcastOperators => BroadcastOperators_POP_FEOTS 
      PROCEDURE :: ScatterSolution    => ScatterSolution_POP_FEOTS 
      PROCEDURE :: GatherSolution     => GatherSolution_POP_FEOTS 
      PROCEDURE :: ScatterMask        => ScatterMask_POP_FEOTS 
      PROCEDURE :: ScatterSource      => ScatterSource_POP_FEOTS 
#endif

      PROCEDURE :: MapAllToDOF      => MapAllToDOF_POP_FEOTS
      PROCEDURE :: MapTracerToDOF   => MapTracerToDOF_POP_FEOTS
      PROCEDURE :: MapHSetToDOF     => MapHSetToDOF_POP_FEOTS
      PROCEDURE :: MapSourceToDOF   => MapSourceToDOF_POP_FEOTS
      PROCEDURE :: MapTracerFromDOF => MapTracerFromDOF_POP_FEOTS

      PROCEDURE :: ForwardStep         => ForwardStep_POP_FEOTS
      PROCEDURE :: ForwardStepEuler    => ForwardStepEuler_POP_FEOTS
      PROCEDURE :: ForwardStepAB2      => ForwardStepAB2_POP_FEOTS
      PROCEDURE :: ForwardStepAB3      => ForwardStepAB3_POP_FEOTS
      PROCEDURE :: CycleIntegrationAB3 => CycleIntegrationAB3_POP_FEOTS
      PROCEDURE :: DotProduct          => DotProduct_POP_FEOTS
      PROCEDURE :: JacobianAction      => JacobianAction_POP_FEOTS
      PROCEDURE :: SolveGMRES       => SolveGMRES_POP_FEOTS
      PROCEDURE :: JFNK             => JFNK_POP_FEOTS
!      PROCEDURE :: FillDiagnostics  => FillDiagnostics_POP_FEOTS

!      PROCEDURE :: ReadSparseConnectivity  => ReadSparseConnectivity_POP_FEOTS
!      PROCEDURE :: ReadMatrixData          => ReadMatrixData_POP_FEOTS
!      PROCEDURE :: WriteTecplot            => WriteTecplot_POP_FEOTS
!      PROCEDURE :: WritePickup             => WritePickup_POP_FEOTS
!      PROCEDURE :: ReadPickup              => ReadPickup_POP_FEOTS
!      PROCEDURE :: WriteDiagnosticsTecplot => WriteDiagnosticsTecplot_POP_FEOTS

   END TYPE POP_FEOTS

#ifdef HAVE_OPENMP
! These arrays are used by the Adams-Bashforth integrators. In order to reduce
! overhead associated with OpenMP and maintain correctness, these arrays need
! to be shared amongst the OpenMP threads. When declared locally within the 
! subroutine, as in the serial code, these arrays become thread private and each
! thread only sees partial updates of each array. When declared here, these
! arrays have global scope and become shared amongst OpenMP threads.
   REAL(prec), ALLOCATABLE :: trm1(:,:)
   REAL(prec), ALLOCATABLE :: trm2(:,:)
   REAL(prec), ALLOCATABLE :: vol(:)
   REAL(prec), ALLOCATABLE :: weightedTracers(:,:)
   REAL(prec), ALLOCATABLE :: dCdt(:,:)
   REAL(prec), ALLOCATABLE :: dVdt(:)

   ! Function JacobianAction
   REAL(prec)              :: scFac, vmag, e
#endif


CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE Build_POP_FEOTS( this, myRank, nProcs ) 
 ! S/R Build
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( POP_FEOTS ), INTENT(out) :: this
   INTEGER, INTENT(in)             :: myRank, nProcs
   ! Local
   INTEGER    :: nX, nY, nZ, nTracers, nRow, nPeriods
   INTEGER    :: nDOF, iLayer, iMask, iTracer
   REAL(prec) :: opPeriod, dt
   INTEGER    :: nOps, i, j, k, m, stencilSize, fUnit
   INTEGER    :: trackingVar
   INTEGER, ALLOCATABLE :: nEl(:)
   CHARACTER(200) :: oceanStateFile
   CHARACTER(20)  :: trackingVar_Char
   

      PRINT*, 'S/R : Build_POP_FEOTS : Start...'

      CALL this % params  % Build( )

      IF( this % params % StencilType == LaxWendroff )THEN
         stencilSize = 7
      ELSE
         STOP 'Bad Stencil!'
      ENDIF

      IF( this % params % Regional )THEN
         IF( TRIM(this % params % maskfile) == '' )THEN
            CALL this % mesh % Load( TRIM(this % params % RegionalMeshFile) )
            CALL this % regionalMaps % ReadPickup( TRIM(this % params % regionalOperatorDirectory)//'mappings', maskProvided=.FALSE. )
         ELSE
            CALL this % mesh % LoadWithMask( TRIM(this % params % RegionalMeshFile) )
            CALL this % regionalMaps % ReadPickup( TRIM(this % params % regionalOperatorDirectory)//'mappings', maskProvided=.TRUE. )
         ENDIF

         this % mesh % DOFtoIJK = this % regionalMaps % dofToLocalIJK
         DO m = 1, this % regionalMaps % nCells
            i = this % regionalMaps % dofToLocalIJK(1,m)
            j = this % regionalMaps % dofToLocalIJK(2,m)
            k = this % regionalMaps % dofToLocalIJK(3,m)
            this % mesh % IJKtoDOF(i,j,k) = m
            this % mesh % DOFtoIJK(1,m)   = i
            this % mesh % DOFtoIJK(2,m)   = j
            this % mesh % DOFtoIJK(3,m)   = k
         ENDDO
         nDOF = this % regionalMaps % nCells
         this % mesh % nDOF = nDOF
      ELSE
         CALL this % mesh % Load( TRIM(this % params % meshFile) )
         nDOF = this % mesh % nDOF
      ENDIF

      IF( this % params % TracerModel == DyeModel )THEN
         nOps   = 2
         ALLOCATE( nEl(1:2) )
         nEl(1) = nDOF*stencilSize !  number of nonzero entries in sparse advection operator
         nEl(2) = nDOF*3 ! number of nonzero entries in sparse vertical diffusion operator
      ELSEIF( this % params % TracerModel == RadionuclideModel )THEN
         nOps = 4
         ALLOCATE( nEl(1:4) )
         nEl(1) = nDOF*stencilSize !  number of nonzero entries in sparse advection operator
         nEl(2) = nDOF*3 ! number of nonzero entries in sparse vertical diffusion operator
         nEl(3) = nDOF*2 ! number of nonzero entries in sparse settling operator
         nEl(4) = nDOF*2 ! number of nonzero entries in sparse scavenging operator
      ELSEIF( this % params % TracerModel == SettlingModel )THEN
         nOps = 3
         ALLOCATE( nEl(1:3) )
         nEl(1) = nDOF*stencilSize !  number of nonzero entries in sparse advection operator
         nEl(2) = nDOF*3 ! number of nonzero entries in sparse vertical diffusion operator
         nEl(3) = nDOF*2 ! number of nonzero entries in sparse settling operator
      ENDIF

      IF( myRank /= 0 )THEN
         this % params % nTracers = 1
      ENDIF
      IF( this % params % WaterMassTagging )THEN
         
         PRINT*, '   Enabling water mass tagging.'
         ! Overwrite the number of tracers
         IF( myRank == 0 )THEN
            this % params % nTracers = this % params % nLayers*this % regionalMaps % nMasks
#ifdef HAVE_MPI
            IF( this % params % nTracers /= nProcs-1 )THEN
               PRINT *, 'Number of tracers (plus one) does not match the number of MPI Ranks.'
               PRINT *, 'nTracers = ', this % params % nTracers
               PRINT *, 'nRanks   = ', nProcs
               STOP 'Stopping!'
            ENDIF
#endif
         ENDIF

         ALLOCATE( this % stateLowerBound(1:this % params % nLayers), &
                   this % stateUpperBound(1:this % params % nLayers) )
         this % stateLowerBound = 0.0_prec
         this % stateUpperBound = 0.0_prec
         
         !IF( LoadForDriver )THEN
            OPEN( UNIT   = NewUnit( fUnit ), &
                  FILE   = 'watermass.config', &
                  FORM   = 'FORMATTED', &
                  STATUS = 'OLD', &
                  ACTION = 'READ' )

            READ( fUnit, '(A20)' ) trackingVar_Char
            PRINT*, '   Tagging water masses with '//TRIM( trackingVar_Char )
            trackingVar = GetFlagforChar( TRIM(trackingVar_Char) ) 
            IF( trackingVar == Temperature ) THEN
               this % stateMask(1) = 1.0_prec
               this % stateMask(2) = 0.0_prec
               this % stateMask(3) = 0.0_prec
            ELSEIF( trackingVar == Salinity ) THEN
               this % stateMask(1) = 0.0_prec
               this % stateMask(2) = 1.0_prec
               this % stateMask(3) = 0.0_prec
            ELSEIF( trackingVar == Density ) THEN
               this % stateMask(1) = 0.0_prec
               this % stateMask(2) = 0.0_prec
               this % stateMask(3) = 1.0_prec
            ELSE
               PRINT*, 'Bad Water Mass Tracking variable.'
               STOP 'STOPPING!'
            ENDIF

            PRINT*, '         Layer |        Lower Bound           |    Upper Bound' 
            DO i = 1, this % params % nLayers
               READ( fUnit, * ) this % stateLowerBound(i), this % stateUpperBound(i)
               PRINT*, i, '   |', this % stateLowerBound(i),'   |', this % stateUpperBound(i)
            ENDDO

            CLOSE( fUnit )
         !ENDIF
      ENDIF
#ifdef HAVE_MPI
      IF( this % params % TracerModel /= DyeModel )THEN
         PRINT*, 'MPI currently only configured for Passive Dye Model.'
         PRINT*, 'Contact Joe at schoonover.numerics@gmail.com for help if needed.'
         STOP 'Stopping!'
      ENDIF
      IF( myRank == 0 )THEN
         IF( this % params % nTracers /= nProcs-1 )THEN
            PRINT *, 'Number of tracers does not match the number of MPI Ranks.'
            PRINT *, 'nTracers = ', this % params % nTracers
            PRINT *, 'nRanks   = ', nProcs
            STOP 'Stopping!'
         ENDIF
      ENDIF
#endif

      ! Allocates space for the solution storage as a 1-D array and allocates
      ! space for the transport operators 
      CALL this % solution % Build( nDOF, nOps, nEl, &
                                    this % params % nOperatorsPerCycle, &
                                    this % params % nTracers, & 
                                    this % params % operatorPeriod, &
                                    this % params % dt )
#ifdef HAVE_MPI
      this % solution % nTracers = 1
#endif

      IF( this % params % Regional ) THEN
         CALL this % solution % LoadSparseConnectivities( &
                                TRIM( this % params % regionalOperatorDirectory)//TRIM(this % params % operatorBasename) )
      ELSE
         CALL this % solution % LoadSparseConnectivities( &
                                TRIM( this % params % feotsOperatorDirectory)//TRIM(this % params % operatorBasename) )
      ENDIF
      IF( this % params % TracerModel == RadionuclideModel .OR. &
          this % params % TracerModel == SettlingModel )THEN 
         PRINT*, '  Setting up Settling Operator'
         CALL this % SetupSettlingOperator( )
      ENDIF

      ! Allocates space for the solution storage on the native mesh 
      CALL this % nativeSol % Build( this % mesh, this % params % nTracers )
      this % nativeSol % mask = 1.0_prec
      IF( this % params % Regional )THEN

         IF( myRank == 0 )THEN
   
            DO iMask = 1, this % regionalMaps % nMasks
               DO iLayer = 1, this % params % nLayers
                     
                  iTracer = iLayer + (iMask-1)*( this % params % nLayers )
                  DO m = 1, this % regionalMaps % bMap(iMask) % nBCells
                     i = this % regionalMaps % dofToLocalIJK(1,this % regionalMaps % bMap(iMask) % boundaryCells(m))
                     j = this % regionalMaps % dofToLocalIJK(2,this % regionalMaps % bMap(iMask) % boundaryCells(m))
                     k = this % regionalMaps % dofToLocalIJK(3,this % regionalMaps % bMap(iMask) % boundaryCells(m))
                     this % nativeSol % mask(i,j,k,iTracer) = 0.0_prec
                  ENDDO
               ENDDO
            ENDDO

         ELSE
            iMask   = (myRank-1)/(this % params % nLayers)+1
            iLayer  = myRank - (iMask-1)*this % params % nLayers
            iTracer = 1
            DO m = 1, this % regionalMaps % bMap(iMask) % nBCells
               i = this % regionalMaps % dofToLocalIJK(1,this % regionalMaps % bMap(iMask) % boundaryCells(m))
               j = this % regionalMaps % dofToLocalIJK(2,this % regionalMaps % bMap(iMask) % boundaryCells(m))
               k = this % regionalMaps % dofToLocalIJK(3,this % regionalMaps % bMap(iMask) % boundaryCells(m))
               this % nativeSol % mask(i,j,k,iTracer) = 0.0_prec
            ENDDO

         ENDIF

      ENDIF

#ifdef HAVE_OPENMP
   ALLOCATE( trm1(1:this % solution % nDOF, 1:this % solution % nTracers), &
             trm2(1:this % solution % nDOF, 1:this % solution % nTracers), &
             vol(1:this % solution % nDOF), &
             weightedTracers(1:this % solution % nDOF, 1:this % solution % nTracers), &
             dCdt(1:this % solution % nDOF, 1:this % solution % nTracers), &
             dVdt(1:this % solution % nDOF) )
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
      CALL this % solution % Trash( )
      IF( this % params % Regional )THEN
         CALL this % regionalMaps % Trash( )
      ENDIF
      CALL this % NativeSol % Trash( )     

#ifdef HAVE_OPENMP
  DEALLOCATE( trm1, &
              trm2, &
              vol, &
              weightedTracers, &
              dCdt, &
              dVdt )
#endif

 END SUBROUTINE Trash_POP_FEOTS
!
 SUBROUTINE SetupSettlingOperator_POP_FEOTS( this )
   IMPLICIT NONE
   CLASS( POP_FEOTS ), INTENT(inout) :: this
   ! Local
   REAL(prec) :: ws, fac
   INTEGER    :: iel, row, col, i, j, k

      ws = this % params % settlingVelocity

      iel = 0
      DO row = 1, this % mesh % nDOF

         i = this % mesh % DOFtoIJK(1,row) ! 
         j = this % mesh % DOFtoIJK(2,row) ! 
         k = this % mesh % DOFtoIJK(3,row) ! vertical level
         
         IF( k > 1 )THEN

            fac = ws/( this % mesh % z(k) - this % mesh % z(k-1) )
          
            col = this % mesh % IJKtoDOF(i,j,k-1)
            iel = iel + 1
            this % solution % transportOps(3) % rowBounds(1,row) = iel
            this % solution % transportOps(3) % A(iel)           = fac              
            this % solution % transportOps(3) % col(iel)         = col ! sub-diagonal   
 
            iel = iel + 1
            this % solution % transportOps(3) % rowBounds(2,row) = iel
            this % solution % transportOps(3) % A(iel)           = -fac              
            this % solution % transportOps(3) % col(iel)         = row ! diagonal              
 
         ELSE

            ! Here, k=1, ie we're at the top-most layer. The flux through the
            ! face is zero
            fac = ws/( this % mesh % z(k)  )
            iel = iel + 1
            this % solution % transportOps(3) % rowBounds(1,row) = iel
            this % solution % transportOps(3) % rowBounds(2,row) = iel
            this % solution % transportOps(3) % A(iel)           = -fac              
            this % solution % transportOps(3) % col(iel)         = row ! diagonal              
 
         ENDIF

      ENDDO

 END SUBROUTINE SetupSettlingOperator_POP_FEOTS
!
!
!==================================================================================================!
!----------------------------------------- Type Specific ------------------------------------------!
!==================================================================================================!
!
!
#ifdef HAVE_MPI
 SUBROUTINE BroadCastOperators_POP_FEOTS( this, myRank, nProcs )
   IMPLICIT NONE
   CLASS( POP_FEOTS ), INTENT(inout) :: this
   INTEGER, INTENT(in)               :: myRank, nProcs
   ! Local
   INTEGER :: iOp

      DO iOp = 1, this % solution % nOps
  
         CALL MPI_Bcast( this % solution % transportOps(iOp) % A, &
                         this % solution % transportOps(iOp) % nElems, &
                         MPI_DOUBLE, &
                         0, MPI_COMM_WORLD )
      ENDDO
 
 END SUBROUTINE BroadcastOperators_POP_FEOTS
!
 SUBROUTINE ScatterSolution_POP_FEOTS( this, myRank, nProcs )
   IMPLICIT NONE
   CLASS( POP_FEOTS ), INTENT(inout) :: this
   INTEGER, INTENT(in)               :: myRank, nProcs
   ! Local
   INTEGER :: mpiErr, i, theStat, recvReq
   INTEGER :: sendReq(1:nProcs-1)
   INTEGER :: theStats(MPI_STATUS_SIZE,1:nProcs-1)


   IF( myRank == 0 )THEN

      DO i = 1, this % params % nTracers
         CALL MPI_ISEND( this % nativeSol % tracer(:,:,:,i), &
                         this % mesh % nX*this % mesh % nY*this % mesh % nZ, &
                         MPI_DOUBLE, &
                         i, 0, MPI_COMM_WORLD, &
                         sendReq(i), mpiErr )
      ENDDO
      CALL MPI_WAITALL( nProcs-1, sendReq, theStats, mpiErr )
   ELSE

      CALL MPI_IRECV( this % nativeSol % tracer(:,:,:,1), &
                     this % mesh % nX*this % mesh % nY*this % mesh % nZ, &
                     MPI_DOUBLE, &
                     0, 0, MPI_COMM_WORLD, &
                     recvReq, mpiErr )
      CALL MPI_WAIT( recvReq, theStat, mpiErr )
   ENDIF
 
 
 END SUBROUTINE ScatterSolution_POP_FEOTS
!
 SUBROUTINE GatherSolution_POP_FEOTS( this, myRank, nProcs )
   IMPLICIT NONE
   CLASS( POP_FEOTS ), INTENT(inout) :: this
   INTEGER, INTENT(in)               :: myRank, nProcs
   ! Local
   INTEGER :: mpiErr, i, theStat, sendReq
   INTEGER :: recvReq(1:nProcs-1)
   INTEGER :: theStats(MPI_STATUS_SIZE,1:nProcs-1)


   !CALL MPI_BARRIER( MPI_COMM_WORLD, mpiErr )
      
   IF( myRank == 0 )THEN

      DO i = 1, this % params % nTracers
         CALL MPI_IRECV( this % nativeSol % tracer(:,:,:,i), &
                         this % mesh % nX*this % mesh % nY*this % mesh % nZ, &
                         MPI_DOUBLE, &
                         i, 0, MPI_COMM_WORLD, &
                         recvReq(i), mpiErr )
      ENDDO
      CALL MPI_WAITALL( nProcs-1, recvReq, theStats, mpiErr )
   ELSE

      CALL MPI_ISEND( this % nativeSol % tracer(:,:,:,1), &
                     this % mesh % nX*this % mesh % nY*this % mesh % nZ, &
                     MPI_DOUBLE, &
                     0, 0, MPI_COMM_WORLD, &
                     sendReq, mpiErr )
      CALL MPI_WAIT( sendReq, theStat, mpiErr )
   ENDIF
 
 END SUBROUTINE GatherSolution_POP_FEOTS
!
 SUBROUTINE ScatterMask_POP_FEOTS( this, myRank, nProcs )
   IMPLICIT NONE
   CLASS( POP_FEOTS ), INTENT(inout) :: this
   INTEGER, INTENT(in)               :: myRank, nProcs
   ! Local
   INTEGER :: mpiErr, i, theStat, recvReq
   INTEGER :: sendReq(1:nProcs-1)
   INTEGER :: theStats(MPI_STATUS_SIZE,1:nProcs-1)


   IF( myRank == 0 )THEN

      PRINT*, myRank, this % params % nTracers
      DO i = 1, this % params % nTracers
         CALL MPI_ISEND( this % solution % mask(:,i), &
                         this % solution % nDOF, &
                         MPI_DOUBLE, &
                         i, 0, MPI_COMM_WORLD, &
                         sendReq(i), mpiErr )
      ENDDO
      CALL MPI_WAITALL( nProcs-1, sendReq, theStats, mpiErr )
   ELSE

      CALL MPI_IRECV( this % solution % mask(:,1), &
                      this % solution % nDOF, &
                      MPI_DOUBLE, &
                      0, 0, MPI_COMM_WORLD, &
                      recvReq, mpiErr )
      CALL MPI_WAIT( recvReq, theStat, mpiErr )
   ENDIF

 END SUBROUTINE ScatterMask_POP_FEOTS
!
 SUBROUTINE ScatterSource_POP_FEOTS( this, myRank, nProcs )
   IMPLICIT NONE
   CLASS( POP_FEOTS ), INTENT(inout) :: this
   INTEGER, INTENT(in)               :: myRank, nProcs
   ! Local
   INTEGER :: mpiErr, i, theStat, recvReq
   INTEGER :: sendReq(1:nProcs-1)
   INTEGER :: theStats(MPI_STATUS_SIZE,1:nProcs-1)


   IF( myRank == 0 )THEN

      PRINT*, myRank, this % params % nTracers
      DO i = 1, this % params % nTracers
         CALL MPI_ISEND( this % solution % source(:,i), &
                         this % solution % nDOF, &
                         MPI_DOUBLE, &
                         i, 0, MPI_COMM_WORLD, &
                         sendReq(i), mpiErr )
      ENDDO
      CALL MPI_WAITALL( nProcs-1, sendReq, theStats, mpiErr )
   ELSE

      CALL MPI_IRECV( this % solution % source(:,1), &
                      this % solution % nDOF, &
                      MPI_DOUBLE, &
                      0, 0, MPI_COMM_WORLD, &
                      recvReq, mpiErr )
      CALL MPI_WAIT( recvReq, theStat, mpiErr )
   ENDIF
   IF( myRank == 0 )THEN

      PRINT*, myRank, this % params % nTracers
      DO i = 1, this % params % nTracers
         CALL MPI_ISEND( this % solution % rFac(:,i), &
                         this % solution % nDOF, &
                         MPI_DOUBLE, &
                         i, 0, MPI_COMM_WORLD, &
                         sendReq(i), mpiErr )
      ENDDO
      CALL MPI_WAITALL( nProcs-1, sendReq, theStats, mpiErr )
   ELSE

      CALL MPI_IRECV( this % solution % rFac(:,1), &
                      this % solution % nDOF, &
                      MPI_DOUBLE, &
                      0, 0, MPI_COMM_WORLD, &
                      recvReq, mpiErr )
      CALL MPI_WAIT( recvReq, theStat, mpiErr )
   ENDIF

 END SUBROUTINE ScatterSource_POP_FEOTS
#endif
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

      DO i = 1, this % params % nTracers 
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

      DO i = 1, this % params % nTracers
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

      DO i = 1, this % params % nTracers
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
 
      DO i = 1, this % params % nTracers
         this % nativeSol % tracer(:,:,:,i) = this % mesh % MapFromDOFtoIJK( this % solution % tracers(:,i) )
      ENDDO
      this % nativeSol % volume = this % mesh % MapFromDOFtoIJK( this % solution % volume )
   
 END SUBROUTINE MapTracerFromDOF_POP_FEOTS
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
#ifndef HAVE_OPENMP
   REAL(prec) :: trm1(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: trm2(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: vol(1:this % solution % nDOF)
   REAL(prec) :: weightedTracers(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: dCdt(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: dVdt(1:this % solution % nDOF)
#endif
   INTEGER         :: i, j, k, m, iT, dof, iTracer, iLayer, iMask
   LOGICAL         :: operatorsSwapped
   CHARACTER(400)  :: oceanStateFile
   CHARACTER(5)    :: fileIDChar
   REAL(prec)      :: trackingVar

      DO iT = 1, nTimeSteps

         !$OMP MASTER
         IF( this % params % Regional ) THEN
            CALL this % solution % CheckForNewOperator( tn, &
                                TRIM( this % params % regionalOperatorDirectory)//TRIM(this % params % operatorBasename), &
                                operatorsSwapped )
         ELSE
            CALL this % solution % CheckForNewOperator( tn, &
                                TRIM( this % params % feotsOperatorDirectory)//TRIM(this % params % operatorBasename), &
                                operatorsSwapped )
         ENDIF

         IF( this % params % waterMassTagging .AND. operatorsSwapped)THEN

            WRITE(fileIDChar, '(I5.5)' ) this % solution % currentPeriod + 1 
            oceanStateFile = TRIM(this % params % regionalOperatorDirectory)//'Ocean.'//fileIDChar//'.nc'
            PRINT*, 'Loading new ocean state : ', TRIM(oceanStateFile)
            CALL this % nativeSol % LoadOceanState( this % mesh, TRIM(oceanStateFile) )

            ! Set any prescribed cells here
#ifndef HAVE_MPI
            DO iMask = 1, this % regionalMaps % nMasks
               DO iLayer = 1, this % params % nLayers
                  
                  iTracer = iLayer + (iMask-1)*( this % params % nLayers )
#else
                  iMask   = (myRank-1)/(this % params % nLayers)+1
                  iLayer  = myRank  - (iMask-1)*this % params % nLayers
                  iTracer = 1
#endif

                  DO m = 1, this % regionalMaps % bMap(iMask) % nPCells

                     dof = this % regionalMaps % bMap(iMask) % prescribedCells(m)
                     i   = this % regionalMaps % dofToLocalIJK(1,dof)
                     j   = this % regionalMaps % dofToLocalIJK(2,dof)
                     k   = this % regionalMaps % dofToLocalIJK(3,dof)

                     trackingVar = this % nativeSol % temperature(i,j,k)*this % statemask(1) +&
                                   this % nativeSol % salinity(i,j,k)*this % statemask(2) +&
                                   this % nativeSol % density(i,j,k)*this % statemask(3)                 
                  
                     this % solution % tracers(dof,iTracer) = 0.0_prec ! Reset prescribed values

                     IF( trackingVar >= this % stateLowerBound(iLayer) .AND. &
                         trackingVar < this % stateUpperBound(iLayer) )THEN
                        this % solution % tracers(dof,iTracer) = 1.0_prec
                     ENDIF

                  ENDDO
#ifndef HAVE_MPI
               ENDDO 
            ENDDO
#endif
         ENDIF
         !$OMP END MASTER

         !$OMP DO COLLAPSE(2)
         DO m = 1, this % solution % nTracers
            DO i = 1, this % solution % nDOF
               weightedTracers(i,m) = this % solution % tracers(i,m)
               trm1(i,m) = this % solution % tracers(i,m)
            ENDDO
         ENDDO
         !$OMP ENDDO


         !$OMP BARRIER
         CALL this % solution % CalculateTendency( weightedTracers, tn, this % params % TracerModel, dCdt, dVdt )
         !$OMP BARRIER

         ! Forward Step the volume 

#ifdef VOLUME_CORRECTION
         !$OMP DO
         DO i = 1, this % solution % nDOF
            vol(i) = this % solution % volume(i) + this % params % dt*dVdt(i)
         ENDDO
         !$OMP ENDDO
#endif

         ! Forward step the tracers with the volume correction
         !$OMP DO COLLAPSE(2)
         DO m = 1, this % solution % nTracers
            DO i = 1, this % solution % nDOF
         !   tracers(:,m)  = (1.0_prec/(1.0_prec+vol))*( (1.0_prec + this % solution % volume )*tracers(:,m) + dt*dCdt(:,m) )
               this % solution % tracers(i,m)  = this % solution % tracers(i,m) + this % params % dt*dCdt(i,m)
            ENDDO
         ENDDO
         !$OMP ENDDO

         ! Store the volume
#ifdef VOLUME_CORRECTION
         !$OMP DO
         DO i = 1, this % solution % nDOF
            this % solution % volume(i) = vol(i)
         ENDDO
         !$OMP ENDDO
#endif

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
#ifndef HAVE_OPENMP
   REAL(prec) :: trm1(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: trm2(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: vol(1:this % solution % nDOF)
   REAL(prec) :: weightedTracers(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: dCdt(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: dVdt(1:this % solution % nDOF)
#endif
   INTEGER         :: i, j, k, m, iT, dof, iTracer, iLayer, iMask
   LOGICAL         :: operatorsSwapped
   CHARACTER(400)  :: oceanStateFile
   CHARACTER(5)    :: fileIDChar
   REAL(prec)      :: trackingVar

      DO iT = 1, nTimeSteps

         !$OMP MASTER
         IF( this % params % Regional ) THEN
            CALL this % solution % CheckForNewOperator( tn, &
                                TRIM( this % params % regionalOperatorDirectory)//TRIM(this % params % operatorBasename), &
                                operatorsSwapped )
         ELSE
            CALL this % solution % CheckForNewOperator( tn, &
                                TRIM( this % params % feotsOperatorDirectory)//TRIM(this % params % operatorBasename), &
                                operatorsSwapped )
         ENDIF

         IF( this % params % waterMassTagging .AND. operatorsSwapped)THEN

            WRITE(fileIDChar, '(I5.5)' ) this % solution % currentPeriod + 1 
            oceanStateFile = TRIM(this % params % regionalOperatorDirectory)//'Ocean.'//fileIDChar//'.nc'
            PRINT*, 'Loading new ocean state : ', TRIM(oceanStateFile)
            CALL this % nativeSol % LoadOceanState( this % mesh, TRIM(oceanStateFile) )

            ! Set any prescribed cells here
#ifndef HAVE_MPI
            DO iMask = 1, this % regionalMaps % nMasks
               DO iLayer = 1, this % params % nLayers
                  
                  iTracer = iLayer + (iMask-1)*( this % params % nLayers )
#else
                  iMask   = (myRank-1)/(this % params % nLayers)+1
                  iLayer  = myRank  - (iMask-1)*this % params % nLayers
                  iTracer = 1
#endif

                  DO m = 1, this % regionalMaps % bMap(iMask) % nPCells

                     dof = this % regionalMaps % bMap(iMask) % prescribedCells(m)
                     i   = this % regionalMaps % dofToLocalIJK(1,dof)
                     j   = this % regionalMaps % dofToLocalIJK(2,dof)
                     k   = this % regionalMaps % dofToLocalIJK(3,dof)

                     trackingVar = this % nativeSol % temperature(i,j,k)*this % statemask(1) +&
                                   this % nativeSol % salinity(i,j,k)*this % statemask(2) +&
                                   this % nativeSol % density(i,j,k)*this % statemask(3)                 
                  
                     this % solution % tracers(dof,iTracer) = 0.0_prec ! Reset prescribed values

                     IF( trackingVar >= this % stateLowerBound(iLayer) .AND. &
                         trackingVar < this % stateUpperBound(iLayer) )THEN
                        this % solution % tracers(dof,iTracer) = 1.0_prec
                     ENDIF

                  ENDDO
#ifndef HAVE_MPI
               ENDDO 
            ENDDO
#endif
         ENDIF
         !$OMP END MASTER

         SELECT CASE (iT)

            CASE(1) ! First order Euler
               !$OMP DO COLLAPSE(2)
               DO m = 1, this % solution % nTracers
                  DO i = 1, this % solution % nDOF
                     weightedTracers(i,m) = this % solution % tracers(i,m)
                     trm1(i,m) = this % solution % tracers(i,m)
                  ENDDO
               ENDDO
               !$OMP ENDDO

            CASE DEFAULT! Second Order Adams Bashforth
               !$OMP DO COLLAPSE(2)
               DO m = 1, this % solution % nTracers
                  DO i = 1, this % solution % nDOF
                     weightedTracers(i,m) = (3.0_prec*this % solution % tracers(i,m) - trm1(i,m))*0.5_prec
                     trm2(i,m) = trm1(i,m)
                     trm1(i,m) = this % solution % tracers(i,m)
                  ENDDO
               ENDDO
               !$OMP ENDDO

         END SELECT

         !$OMP BARRIER
         CALL this % solution % CalculateTendency( weightedTracers, tn, this % params % TracerModel, dCdt, dVdt )
         !$OMP BARRIER

         ! Forward Step the volume 

#ifdef VOLUME_CORRECTION
         !$OMP DO
         DO i = 1, this % solution % nDOF
            vol(i) = this % solution % volume(i) + this % params % dt*dVdt(i)
         ENDDO
         !$OMP ENDDO
#endif

         ! Forward step the tracers with the volume correction
         !$OMP DO COLLAPSE(2)
         DO m = 1, this % solution % nTracers
            DO i = 1, this % solution % nDOF
         !   tracers(:,m)  = (1.0_prec/(1.0_prec+vol))*( (1.0_prec + this % solution % volume )*tracers(:,m) + dt*dCdt(:,m) )
               this % solution % tracers(i,m)  = this % solution % tracers(i,m) + this % params % dt*dCdt(i,m)
            ENDDO
         ENDDO
         !$OMP ENDDO

         ! Store the volume
#ifdef VOLUME_CORRECTION
         !$OMP DO
         DO i = 1, this % solution % nDOF
            this % solution % volume(i) = vol(i)
         ENDDO
         !$OMP ENDDO
#endif

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
#ifndef HAVE_OPENMP
   REAL(prec) :: trm1(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: trm2(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: vol(1:this % solution % nDOF)
   REAL(prec) :: weightedTracers(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: dCdt(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: dVdt(1:this % solution % nDOF)
#endif
   INTEGER         :: i, j, k, m, iT, dof, iTracer, iLayer, iMask
   LOGICAL         :: operatorsSwapped
   CHARACTER(400)  :: oceanStateFile
   CHARACTER(5)    :: fileIDChar
   REAL(prec)      :: trackingVar

      DO iT = 1, nTimeSteps

         !$OMP MASTER
         IF( this % params % Regional ) THEN
            CALL this % solution % CheckForNewOperator( tn, &
                                TRIM( this % params % regionalOperatorDirectory)//TRIM(this % params % operatorBasename), &
                                operatorsSwapped )
         ELSE
            CALL this % solution % CheckForNewOperator( tn, &
                                TRIM( this % params % feotsOperatorDirectory)//TRIM(this % params % operatorBasename), &
                                operatorsSwapped )
         ENDIF

         IF( this % params % waterMassTagging .AND. operatorsSwapped)THEN

            WRITE(fileIDChar, '(I5.5)' ) this % solution % currentPeriod + 1 
            oceanStateFile = TRIM(this % params % regionalOperatorDirectory)//'Ocean.'//fileIDChar//'.nc'
            PRINT*, 'Loading new ocean state : ', TRIM(oceanStateFile)
            CALL this % nativeSol % LoadOceanState( this % mesh, TRIM(oceanStateFile) )

            ! Set any prescribed cells here
#ifndef HAVE_MPI
            DO iMask = 1, this % regionalMaps % nMasks
               DO iLayer = 1, this % params % nLayers
                  
                  iTracer = iLayer + (iMask-1)*( this % params % nLayers )
#else
                  iMask   = (myRank-1)/(this % params % nLayers)+1
                  iLayer  = myRank  - (iMask-1)*this % params % nLayers
                  iTracer = 1
#endif

                  DO m = 1, this % regionalMaps % bMap(iMask) % nPCells

                     dof = this % regionalMaps % bMap(iMask) % prescribedCells(m)
                     i   = this % regionalMaps % dofToLocalIJK(1,dof)
                     j   = this % regionalMaps % dofToLocalIJK(2,dof)
                     k   = this % regionalMaps % dofToLocalIJK(3,dof)

                     trackingVar = this % nativeSol % temperature(i,j,k)*this % statemask(1) +&
                                   this % nativeSol % salinity(i,j,k)*this % statemask(2) +&
                                   this % nativeSol % density(i,j,k)*this % statemask(3)                 
                  
                     this % solution % tracers(dof,iTracer) = 0.0_prec ! Reset prescribed values

                     IF( trackingVar >= this % stateLowerBound(iLayer) .AND. &
                         trackingVar < this % stateUpperBound(iLayer) )THEN
                        this % solution % tracers(dof,iTracer) = 1.0_prec
                     ENDIF

                  ENDDO
#ifndef HAVE_MPI
               ENDDO 
            ENDDO
#endif
         ENDIF
         !$OMP END MASTER

         SELECT CASE (iT)

            CASE(1) ! First order Euler
               !$OMP DO COLLAPSE(2)
               DO m = 1, this % solution % nTracers
                  DO i = 1, this % solution % nDOF
                     weightedTracers(i,m) = this % solution % tracers(i,m)
                     trm1(i,m) = this % solution % tracers(i,m)
                  ENDDO
               ENDDO
               !$OMP ENDDO

            CASE (2) ! Second Order Adams Bashforth
               !$OMP DO COLLAPSE(2)
               DO m = 1, this % solution % nTracers
                  DO i = 1, this % solution % nDOF
                     weightedTracers(i,m) = (3.0_prec*this % solution % tracers(i,m) - trm1(i,m))*0.5_prec
                     trm2(i,m) = trm1(i,m)
                     trm1(i,m) = this % solution % tracers(i,m)
                  ENDDO
               ENDDO
               !$OMP ENDDO

            CASE DEFAULT ! Third Order Adams Bashforth
               !$OMP DO COLLAPSE(2)
               DO m = 1, this % solution % nTracers
                  DO i = 1, this % solution % nDOF
                     weightedTracers(i,m) = (23.0_prec*this % solution % tracers(i,m) - 16.0_prec*trm1(i,m) + 5.0_prec*trm2(i,m))/12.0_prec
                     trm2(i,m) = trm1(i,m)
                     trm1(i,m) = this % solution % tracers(i,m)
                  ENDDO
               ENDDO
               !$OMP ENDDO

         END SELECT

         !$OMP BARRIER
         CALL this % solution % CalculateTendency( weightedTracers, tn, this % params % TracerModel, dCdt, dVdt )
         !$OMP BARRIER

         ! Forward Step the volume 

#ifdef VOLUME_CORRECTION
         !$OMP DO
         DO i = 1, this % solution % nDOF
            vol(i) = this % solution % volume(i) + this % params % dt*dVdt(i)
         ENDDO
         !$OMP ENDDO
#endif

         ! Forward step the tracers with the volume correction
         !$OMP DO COLLAPSE(2)
         DO m = 1, this % solution % nTracers
            DO i = 1, this % solution % nDOF
         !   tracers(:,m)  = (1.0_prec/(1.0_prec+vol))*( (1.0_prec + this % solution % volume )*tracers(:,m) + dt*dCdt(:,m) )
               this % solution % tracers(i,m)  = this % solution % tracers(i,m) + this % params % dt*dCdt(i,m)
            ENDDO
         ENDDO
         !$OMP ENDDO

         ! Store the volume
#ifdef VOLUME_CORRECTION
         !$OMP DO
         DO i = 1, this % solution % nDOF
            this % solution % volume(i) = vol(i)
         ENDDO
         !$OMP ENDDO
#endif

         !$OMP MASTER
         tn = tn + this % params % dt
         !$OMP END MASTER

      ENDDO
   
 END SUBROUTINE ForwardStepAB3_POP_FEOTS
!
 SUBROUTINE CycleIntegrationAB3_POP_FEOTS( this )
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
   ! Local
#ifndef HAVE_OPENMP
   REAL(prec) :: trm1(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: trm2(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: vol(1:this % solution % nDOF)
   REAL(prec) :: weightedTracers(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: dCdt(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: dVdt(1:this % solution % nDOF)
#endif
   REAL(prec) :: tn, dt
   INTEGER    :: nPeriods, nSteps, i, j, k, m 
   CHARACTER(5) :: periodChar
   CHARACTER(500) :: fileBase

      nSteps = INT( this % solution % opPeriod/ this % solution % dt )

      DO j = 1, this % params % nOperatorsPerCycle

         ! Load in the transport operators
         !$OMP MASTER
         WRITE( periodChar, '(I5.5)' ) j
         IF( this % params % Regional )THEN
            fileBase = TRIM( this % params % regionalOperatorDirectory)//TRIM(this % params % operatorBasename)
         ELSE
            fileBase = TRIM( this % params % feotsOperatorDirectory)//TRIM(this % params % operatorBasename)
         ENDIF
         CALL this % solution % transportOps(1) % ReadMatrixData(TRIM(fileBase)//'_advect.'//periodChar )
         CALL this % solution % transportOps(2) % ReadMatrixData(TRIM(fileBase)//'_vdiffu.'//periodChar )
         !$OMP END MASTER 

         DO i = 1, nSteps

            tn = REAL(i-1,prec)*dt

            SELECT CASE (i)

               CASE(1) ! First order Euler
                  !$OMP DO COLLAPSE(2)
                  DO m = 1, this % solution % nTracers
                     DO k = 1, this % solution % nDOF
                        weightedTracers(k,m) = this % solution % tracers(k,m)
                        trm1(k,m) = this % solution % tracers(k,m)
                     ENDDO
                  ENDDO
                  !$OMP ENDDO
               CASE(2) ! Second Order Adams Bashforth
                  !$OMP DO COLLAPSE(2)
                  DO m = 1, this % solution % nTracers
                     DO k = 1, this % solution % nDOF
                        weightedTracers(k,m) = (3.0_prec*this % solution % tracers(k,m) - trm1(k,m))*0.5_prec
                        trm2(k,m) = trm1(k,m)
                        trm1(k,m) = this % solution % tracers(k,m)
                     ENDDO
                  ENDDO
                  !$OMP ENDDO
               CASE DEFAULT ! Third Order Adams Bashforth
                  !$OMP DO COLLAPSE(2)
                  DO m = 1, this % solution % nTracers
                     DO k = 1, this % solution % nDOF
                        weightedTracers(k,m) = (23.0_prec*this % solution % tracers(k,m) - 16.0_prec*trm1(k,m) + 5.0_prec*trm2(k,m))/12.0_prec
                        trm2(k,m) = trm1(k,m)
                        trm1(k,m) = this % solution % tracers(k,m)
                     ENDDO
                  ENDDO
                  !$OMP ENDDO
            END SELECT

            !$OMP BARRIER
            CALL this % solution % CalculateTendency( weightedTracers, tn, this % params % TracerModel, dCdt, dVdt )
            !$OMP BARRIER

            ! Forward Step the volume 
            !$OMP DO
            DO k = 1, this % solution % nDOF 
               vol(k) = this % solution % volume(k) + this % params % dt*dVdt(k)
            ENDDO
            !$OMP ENDDO

            ! Forward step the tracers with the volume correction
            !$OMP DO COLLAPSE(2)
            DO m = 1, this % solution % nTracers
               DO k = 1, this % solution % nDOF
            !   tracers(:,m)  =(1.0_prec/(1.0_prec+vol))*( (1.0_prec + this % solution % volume )*tracers(:,m) + dt*dCdt(:,m) )
                  this % solution % tracers(k,m)  = this % solution % tracers(k,m) + this % params % dt*dCdt(k,m)
               ENDDO
            ENDDO
            !$OMP ENDDO

            ! Store the volume
            !$OMP DO
            DO k = 1, this % solution % nDOF
               this % solution % volume(k) = vol(k) 
            ENDDO
            !$OMP ENDDO

         ENDDO

         ! This last step ensures we end on t=tCycle
         tn = REAL(nSteps,prec)*this % params % dt
         dt = this % params % operatorPeriod - tn
         CALL this % solution % CalculateTendency( this % solution % tracers, tn, this % params % TracerModel, dCdt, dVdt )
         
         !$OMP DO
         DO k = 1, this % solution % nDOF
            vol(k) = this % solution % volume(k) + dt*dVdt(k)
         ENDDO
         !$OMP ENDDO

         !$OMP DO COLLAPSE(2)
         DO m = 1, this % solution % nTracers
            DO k = 1, this % solution % nDOF
           ! tracers(:,m)  =(1.0_prec/(1.0_prec+vol))*( (1.0_prec + this % solution % volume )*tracers(:,m) + dt*dCdt(:,m) )
               this % solution % tracers(k,m)  = this % solution % tracers(k,m) + dt*dCdt(k,m)
            ENDDO
         ENDDO
         !$OMP ENDDO

         !$OMP DO
         DO k = 1, this % solution % nDOF
            this % solution % volume(k) = vol(k)
         ENDDO
         !$OMP ENDDO

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

      xdoty = 0.0_prec
      DO j = 1, this % solution % nTracers
         DO i = 1, this % solution % nDOF
            xdoty = xdoty + x(i,j)*y(i,j)
         ENDDO
      ENDDO

 END FUNCTION DotProduct_POP_FEOTS
!
! ================================================================================================ !
!                                  Equilibration Routines
! ================================================================================================ !
!
 FUNCTION JacobianAction_POP_FEOTS( this, x, Gx, v ) RESULT( Jv )
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
   ! Local
   INTEGER           :: i, j
#ifndef HAVE_OPENMP
   REAL(prec)        :: scFac, vmag, e
#endif

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

      CALL this % CycleIntegrationAB3( )
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
 SUBROUTINE SolveGMRES_POP_FEOTS( this, x, Gx, dx, resi, ioerr )
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
#ifndef HAVE_OPENMP
   REAL(prec) :: e, vmag, scFac
#endif 

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

            CALL this % CycleIntegrationAB3( )
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

            CALL this % CycleIntegrationAB3( )
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
 SUBROUTINE JFNK_POP_FEOTS( this )
 !
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( POP_FEOTS ), INTENT(inout) :: this
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

      CALL this % CycleIntegrationAB3( )
      Gx = this % solution % tracers - x
      Gx0Mag = SQRT( this % DotProduct( Gx, Gx ) )

      DO iter = 1, iterMax
         PRINT*, 'S/R JFNK : Iterate = ', iter

         ! The first call to "CycleIntegration" calculates G(x_i),
         ! which is the RHS of the Linearized fixed point problem
         !    J*dx = -Gx

         
         PRINT*, 'S/R JFNK : Call to SolveGMRES '
         CALL this % SolveGMRES( x, Gx, dx, residual, ioerr )
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
         CALL this % CycleIntegrationAB3( )
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
                                                              'Tracer.pickup.nc' )
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
                                                           'Tracer.equilibrium.nc' )
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
                                                           'Tracer.pickup.nc' )
         CALL this % nativeSol % WriteNetCDFRecord( this % mesh, 1 )
         CALL this % nativeSol % FinalizeNetCDF( )


         PRINT*, ' Module POP_FEOTS_Class.f90 : S/R JFNK_POP_FEOTS : Solution not found in ', iter, ' iterates.'
         PRINT*, ' Final Residual : ', GxMag
         CLOSE(fUnit)
      ENDIF
         

 END SUBROUTINE JFNK_POP_FEOTS
!
END MODULE POP_FEOTS_Class
