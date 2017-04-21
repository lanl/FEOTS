! POP_FEOTS_Class.f90
! 
! Copyright 2016 Joseph Schoonover, Los Alamos National Laboratory (jschoonover@lanl.gov) 
! All rights reserved. 
! 
! POP_FEOTS_Class.f90 is part of the Fast Equilibration of Ocean Tracers Software (FEOTS). 
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

      CONTAINS

      PROCEDURE :: Build => Build_POP_FEOTS
      PROCEDURE :: Trash => Trash_POP_FEOTS

      PROCEDURE :: SetupSettlingOperator => SetupSettlingOperator_POP_FEOTS

      PROCEDURE :: MapAllToDOF      => MapAllToDOF_POP_FEOTS
      PROCEDURE :: MapTracerToDOF   => MapTracerToDOF_POP_FEOTS
      PROCEDURE :: MapHSetToDOF     => MapHSetToDOF_POP_FEOTS
      PROCEDURE :: MapSourceToDOF   => MapSourceToDOF_POP_FEOTS
      PROCEDURE :: MapTracerFromDOF => MapTracerFromDOF_POP_FEOTS

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


CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE Build_POP_FEOTS( this ) 
 ! S/R Build
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( POP_FEOTS ), INTENT(out) :: this
   ! Local
   INTEGER    :: nX, nY, nZ, nTracers, nRow, nPeriods
   REAL(prec) :: opPeriod, dt
   INTEGER :: nOps, i, j, k, m, stencilSize
   INTEGER, ALLOCATABLE :: nEl(:)

      PRINT*, 'S/R : Build_POP_FEOTS : Start...'

      CALL this % params  % Build( )

      IF( this % params % StencilType == LaxWendroff )THEN
         stencilSize = 7
      ELSE
         STOP 'Bad Stencil!'
      ENDIF

      IF( this % params % Regional )THEN
         CALL this % mesh % Load( TRIM(this % params % RegionalMeshFile) )
         CALL this % regionalMaps % ReadPickup( TRIM(this % params % regionalOperatorDirectory)//'mappings' )

         this % mesh % DOFtoIJK = this % regionalMaps % dofToLocalIJK
         DO m = 1, this % regionalMaps % nCells
            i = this % regionalMaps % dofToLocalIJK(1,m)
            j = this % regionalMaps % dofToLocalIJK(2,m)
            k = this % regionalMaps % dofToLocalIJK(3,m)
            this % mesh % IJKtoDOF(i,j,k) = m
         ENDDO

      ELSE
         CALL this % mesh % Load( TRIM(this % params % meshFile) )
      ENDIF

      IF( this % params % TracerModel == DyeModel )THEN
         nOps   = 2
         ALLOCATE( nEl(1:2) )
         nEl(1) = this % mesh % nDOF*stencilSize !  number of nonzero entries in sparse advection operator
         nEl(2) = this % mesh % nDOF*3 ! number of nonzero entries in sparse vertical diffusion operator
      ELSEIF( this % params % TracerModel == RadionuclideModel )THEN
         nOps = 4
         ALLOCATE( nEl(1:4) )
         nEl(1) = this % mesh % nDOF*stencilSize !  number of nonzero entries in sparse advection operator
         nEl(2) = this % mesh % nDOF*3 ! number of nonzero entries in sparse vertical diffusion operator
         nEl(3) = this % mesh % nDOF*2 ! number of nonzero entries in sparse settling operator
         nEl(4) = this % mesh % nDOF*2 ! number of nonzero entries in sparse scavenging operator
      ELSEIF( this % params % TracerModel == SettlingModel )THEN
         nOps = 3
         ALLOCATE( nEl(1:3) )
         nEl(1) = this % mesh % nDOF*stencilSize !  number of nonzero entries in sparse advection operator
         nEl(2) = this % mesh % nDOF*3 ! number of nonzero entries in sparse vertical diffusion operator
         nEl(3) = this % mesh % nDOF*2 ! number of nonzero entries in sparse settling operator
      ENDIF

      ! Allocates space for the solution storage as a 1-D array and allocates
      ! space for the transport operators 
      CALL this % solution % Build( this % mesh % nDOF, nOps, nEl, &
                                    this % params % nOperatorsPerCycle, &
                                    this % params % nTracers+1, &        ! Always add one tracer for the volume correction
                                    this % params % operatorPeriod, &
                                    this % params % dt )

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
      CALL this % nativeSol % Build( this % mesh, this % params % nTracers+1 )
      this % nativeSol % mask = 1.0_prec
      IF( this % params % Regional )THEN

         DO m = 1, this % regionalMaps % nBCells
            i = this % regionalMaps % dofToLocalIJK(1,this % regionalMaps % boundaryCells(m))
            j = this % regionalMaps % dofToLocalIJK(2,this % regionalMaps % boundaryCells(m))
            k = this % regionalMaps % dofToLocalIJK(3,this % regionalMaps % boundaryCells(m))
            this % nativeSol % mask(i,j,k,:) = 0.0_prec
         ENDDO

      ENDIF

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

      DO i = 1, this % params % nTracers+1 
         this % solution % tracers(:,i) = this % mesh % MapFromIJKtoDOF( this % nativeSol % tracer(:,:,:,i) )
      ENDDO
   
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

      DO i = 1, this % params % nTracers+1
         this % solution % hSetTracers(:,i) = this % mesh % MapFromIJKtoDOF( this % nativeSol % hardSet(:,:,:,i) )
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

      DO i = 1, this % params % nTracers+1
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
 
      DO i = 1, this % params % nTracers+1
         this % nativeSol % tracer(:,:,:,i) = this % mesh % MapFromDOFtoIJK( this % solution % tracers(:,i) )
      ENDDO
   
 END SUBROUTINE MapTracerFromDOF_POP_FEOTS
!
 SUBROUTINE ForwardStepAB3_POP_FEOTS( this, tn, nTimeSteps )
 ! S/R ForwardStepAB3
 !
 !  Identical to "CycleIntegrationAB3", except that the operators are loaded in
 !  as a function of time
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( POP_FEOTS ), INTENT(inout) :: this
   REAL(prec), INTENT(inout)         :: tn
   INTEGER, INTENT(in)               :: nTimeSteps
   ! Local
   REAL(prec) :: trm1(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: trm2(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: vol(1:this % solution % nDOF)
   REAL(prec) :: weightedTracers(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: tracers(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: dCdt(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: dt, opPeriod
   INTEGER    :: nPeriods, nSteps, i, j, k, m 

      dt       = this % solution % dt
      opPeriod = this % solution % opPeriod
      tracers = this % solution % tracers

      DO k = 1, nTimeSteps

         IF( this % params % Regional ) THEN
            CALL this % solution % CheckForNewOperator( tn, &
                                TRIM( this % params % regionalOperatorDirectory)//TRIM(this % params % operatorBasename) )
         ELSE
            CALL this % solution % CheckForNewOperator( tn, &
                                TRIM( this % params % feotsOperatorDirectory)//TRIM(this % params % operatorBasename) )
         ENDIF

         SELECT CASE (k)

            CASE(1) ! First order Euler
               weightedTracers = tracers
               trm1 = tracers

            CASE(2) ! Second Order Adams Bashforth
               weightedTracers = (3.0_prec*tracers - trm1)*0.5_prec
               trm2 = trm1
               trm1 = tracers

            CASE DEFAULT ! Third Order Adams Bashforth
               weightedTracers = (23.0_prec*tracers - 16.0_prec*trm1 + 5.0_prec*trm2)/12.0_prec
               trm2 = trm1
               trm1 = tracers

         END SELECT

         CALL this % solution % CalculateTendency( weightedTracers, tn, this % params % TracerModel, dCdt )

         ! Forward Step the volume 
         vol = tracers(:,this % solution % nTracers) + dt*dCdt(:,this % solution % nTracers)
         ! Forward step the tracers with the volume correction
         DO m = 1, this % solution % nTracers-1
           ! tracers(:,m)  =(1.0_prec/(1.0_prec+vol))*( (1.0_prec + tracers(:,this % solution % nTracers) )*tracers(:,m) + dt*dCdt(:,m) )
            tracers(:,m)  = tracers(:,m) + dt*dCdt(:,m)
         ENDDO
         ! Store the volume
         tracers(:,this % solution % nTracers) = vol
         tn = tn + dt

      ENDDO

      this % solution % tracers = tracers
 
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
   REAL(prec) :: trm1(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: trm2(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: vol(1:this % solution % nDOF)
   REAL(prec) :: weightedTracers(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: tracers(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: dCdt(1:this % solution % nDOF, 1:this % solution % nTracers)
   REAL(prec) :: tn, dt, opPeriod
   INTEGER    :: nPeriods, nSteps, i, j, k, m 
   CHARACTER(5) :: periodChar
   CHARACTER(500) :: fileBase

      dt       = this % solution % dt
      opPeriod = this % solution % opPeriod

      nSteps = INT( this % solution % opPeriod/ this % solution % dt )

      tracers = this % solution % tracers

      DO j = 1, this % params % nOperatorsPerCycle

         ! Load in the transport operators
         WRITE( periodChar, '(I5.5)' ) j
         IF( this % params % Regional )THEN
            fileBase = TRIM( this % params % regionalOperatorDirectory)//TRIM(this % params % operatorBasename)
         ELSE
            fileBase = TRIM( this % params % feotsOperatorDirectory)//TRIM(this % params % operatorBasename)
         ENDIF
         CALL this % solution % transportOps(1) % ReadMatrixData(TRIM(fileBase)//'_advect.'//periodChar )
         CALL this % solution % transportOps(2) % ReadMatrixData(TRIM(fileBase)//'_vdiffu.'//periodChar )
         ! 

         DO i = 1, nSteps

            tn = REAL(i-1,prec)*dt

            SELECT CASE (i)

               CASE(1) ! First order Euler
                  weightedTracers = tracers
                  trm1 = tracers

               CASE(2) ! Second Order Adams Bashforth
                  weightedTracers = (3.0_prec*tracers - trm1)*HALF
                  trm2 = trm1
                  trm1 = tracers

               CASE DEFAULT ! Third Order Adams Bashforth
                  weightedTracers = (23.0_prec*tracers - 16.0_prec*trm1 + 5.0_prec*trm2)/12.0_prec
                  trm2 = trm1
                  trm1 = tracers

            END SELECT

            CALL this % solution % CalculateTendency( weightedTracers, tn, this % params % TracerModel, dCdt )

            ! Forward Step the volume 
            vol = tracers(:,this % solution % nTracers) + dt*dCdt(:,this % solution % nTracers)

            ! Forward step the tracers with the volume correction
            DO m = 1, this % solution % nTracers-1
               tracers(:,m)  =(1.0_prec/(1.0_prec+vol))*( (1.0_prec + tracers(:,this % solution % nTracers) )*tracers(:,m) + dt*dCdt(:,m) )
            ENDDO

            ! Store the volume
            tracers(:,this % solution % nTracers) = vol 

         ENDDO

         ! This last step ensures we end on t=tCycle
         tn = REAL(nSteps,prec)*dt
         dt = opPeriod - tn
         CALL this % solution % CalculateTendency( tracers, tn, this % params % TracerModel, dCdt )
         tracers  = tracers + dt*dCdt
         !this % time = this % time + dt

      ENDDO

      this % solution % tracers = tracers
 
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
   INTEGER           :: i
   REAL(prec)        :: Gxpv(1:this % solution % nDOF,1:this % solution % nTracers)
   REAL(prec)        :: b, e, vmag, scFac
      
      b = this % params % JacobianStepSize
      
      vmag = SQRT(  this % DotProduct( v, v ) )
      DO i = 1, this % solution % nDOF 
         scFac = scFac + b*( abs(x(i,1)) + abs(x(i,2)) )
      ENDDO
      scFac = scFac + b

      e = scFac/( REAL( this % solution % nDOF*this % solution % nTracers, prec)*vmag )

      this % solution % tracers = x + e*v
      CALL this % CycleIntegrationAB3( )
      !      M( x + e*v )         - (x + e*v)
      Gxpv = this % solution % tracers - (x + e*v)

      Jv = (Gxpv - Gx)/e

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
   REAL(prec) :: b, d, g, r0, rc
   

      ioerr = -2
      nIt = this % params % maxItersGMRES
      m   = this % params % mInnerItersGMRES
      TOL = this % params % toleranceGMRES
      
      ! Assume that the initial guess is dx = 0
      r = -Gx
      r0 = sqrt( this % DotProduct( r, r ) )
     
      l = 0 
      resi = ZERO
      resi(l) = r0
      
      dx   = ZERO
      v    = ZERO
      rho  = ZERO
      bhat = ZERO 
      s    = ZERO
      c    = ZERO
      y    = ZERO
   

      DO j = 1,nIt

         b        = sqrt( this % DotProduct( r, r ) )
         v(:,:,1) = r/b
         bhat(1)  = b

         DO i = 1, m
            l = l+1
            nr = i

            ! Applying the precondtioner (right now it is the identity, ie, no preconditioning)
            z(:,:,i) = v(:,:,i)
            !z(:,i) = this % PreconditionInvert( v(:,i) )
            ! The first step in GMRES is to build the orthogonal basis up to order "i"
            ! with the accompanying upper hessenburg matrix.
            w = this % JacobianAction( x, Gx, z(:,:,i) )

            ! The new basis vector is obtained by multiplying the previous basis vector by the matrix
            ! and orthogonalizing wrt to all of the previous basis vectors using a Gram-Schmidt process.
            DO k = 1, i
               h(k,i) = this % DotProduct( v(:,:,k), w )
               w      = w - h(k,i)*v(:,:,k)
            ENDDO

            h(i+1,i) = sqrt( this % DotProduct(w,w) )

            IF( AlmostEqual( h(i+1,i), ZERO )  )THEN
               EXIT
            ENDIF

            v(:,:,i+1) = w/h(i+1,i)
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
            
            IF( rc/r0 <= TOL )THEN
               EXIT
            ENDIF

         ENDDO

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
        
         DO i = 1, this % solution % nTracers
            localz = z(:,i,1:m)
            dx(:,i) = dx(:,i) + MATMUL( localz(:,1:nr), y(1:nr) )
         ENDDO

         IF( rc/r0 <= TOL )THEN
            ioerr = l
            EXIT
         ENDIF
         
         r = -Gx - this % JacobianAction( x, Gx, dx )
         
      ENDDO 

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
      dx = ZERO

      tk = 1

      OPEN( UNIT   = NewUnit(fUnit), &
            FILE   = 'Residual.curve', &
            FORM   = 'FORMATTED', &
            STATUS = 'REPLACE', &
            ACTION = 'WRITE' )
      WRITE( fUnit, * )'#Residual'

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

      ENDDO

      IF( GxMag <= tolerance )THEN
         ! Converged, write the solution !
         PRINT*, ' Module POP_FEOTS_Class.f90 : S/R JFNK_POP_FEOTS : Solution found in ', iter, ' iterates.'
         PRINT*, ' Final Residual : ', GxMag, Gx0Mag
         PRINT*, ' Final Solution Update magnitude : ', dxMag, xMag
         this % solution % tracers = x
         CALL this % MapTracerFromDOF( ) ! Map from DOF to native storage
         CLOSE(fUnit)
      ELSE
         ! Not converged, write a pickup file !
         !WRITE( iterChar, '(I4.4)' ) iterMax + this % params % cycleStart
         this % solution % tracers = x 
         CALL this % MapTracerFromDOF( )
         !


         PRINT*, ' Module POP_FEOTS_Class.f90 : S/R JFNK_POP_FEOTS : Solution not found in ', iter, ' iterates.'
         PRINT*, ' Final Residual : ', GxMag
         CLOSE(fUnit)
      ENDIF
         

 END SUBROUTINE JFNK_POP_FEOTS
!
END MODULE POP_FEOTS_Class
