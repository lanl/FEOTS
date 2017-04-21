! TracerStorage_Class.f90
! 
! Copyright 2016 Joseph Schoonover, Los Alamos National Laboratory (jschoonover@lanl.gov) 
! All rights reserved. 
! 
! TracerStorage_Class.f90 is part of the Fast Equilibration of Ocean Tracers Software (FEOTS). 
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !
 

MODULE TracerStorage_Class


 USE ModelPrecision
 USE ConstantsDictionary
 USE CommonRoutines
 USE CRSMatrix_Class

 IMPLICIT NONE

! ================================ Tracer Storage Description ==================================== !
!
! 15.04.2016
!
! The TracerStorage class stores the transport operators and the tracer fields in the form of 
! sparse matrices and a two dimensional array respectively. This class assumes that multiple
! transport operators are diagnosed to form a "transport cycle". For example, as in Bardin (2014),
! twelve monthly averaged transport operators we diagnosed to give a one year transport cycle.
!
! Attributes and routines are provided for doing a "source-relaxation" for the tracer field and 
! for "hard-setting" tracers.
!
! In addition to the usual OO routines (constructor/destructor and accessors), type specific 
! routines are provided for performing integration of the transport operators over a complete cycle.
!
! It is intended that the user will extend this class to define additional integration schemes and
! any other routines, such as those needed to load in the transport operators.
!
! Attribute description
!   nDOF         -- Integer, the number of degrees of freedom, ie, the number of wet-points in a GCM
!   nPeriods         -- Integer, the number of transport operators that make up a cycle
!   nTracers     -- Integer, the number of tracers that are being worked with
!   opPeriod     -- Real(prec), the length of time over which a single tranport operator is applied;
!                               opPeriod*nPeriods gives the length of time for a transport cycle.
!   dt           -- Real(prec), the time step size for performing integration
!   transportOps -- SparseMatrix, the transport operators stored as a sparse matrices
!   tracers      -- Real(prec), a two dimensional array, indexed (1:nDOF,1:nTracers), the tracers
!   source       -- Real(prec), " ", " ", the source terms
!   rFac         -- Real(prec), " ", " ", the relaxation factor for the source terms
!   mask         -- Real(prec), " ", " ", a mask for hard-setting tracer values
!   hSetTracers  -- Real(prec), " ", " ", the "hard-set" tracer values where the mask is applied
!
! ================================================================================================ !
!

   TYPE TracerStorage 
      INTEGER                        :: nDOF, nPeriods, nOps, nTracers
      INTEGER                        :: currentPeriod
      REAL(prec)                     :: opPeriod, dt
      TYPE( CRSMatrix ), ALLOCATABLE :: transportOps(:) 
      REAL(prec), ALLOCATABLE        :: tracers(:,:)
      REAL(prec), ALLOCATABLE        :: source(:,:), rFac(:,:)
      REAL(prec), ALLOCATABLE        :: mask(:,:), hSetTracers(:,:)

      CONTAINS
      
      PROCEDURE :: Build => Build_TracerStorage
      PROCEDURE :: Trash => Trash_TracerStorage

      PROCEDURE :: LoadSparseConnectivities => LoadSparseConnectivities_TracerStorage
      PROCEDURE :: CheckForNewOperator => CheckForNewOperator_TracerStorage
      PROCEDURE :: MaskField         => MaskField_TracerStorage
      PROCEDURE :: MaskTracers       => MaskTracers_TracerStorage
      PROCEDURE :: CalculateTendency => CalculateTendency_TracerStorage

   END TYPE TracerStorage


!   PRIVATE :: PassiveDyeModel, ParticulateRadioNuclideModel

 CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE Build_TracerStorage( thisStorage, nDOF, nOps, nElems, nPeriods, nTracers, opPeriod, dt )
 ! S/R Build
 !
 !    Allocates memory for the tracer storage class.
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( TracerStorage ), INTENT(out) :: thisStorage
   INTEGER, INTENT(in)                 :: nDOF, nOps
   INTEGER, INTENT(in)                 :: nElems(1:nOps)
   INTEGER, INTENT(in)                 :: nPeriods, nTracers
   REAL(prec), INTENT(in)              :: opPeriod, dt
   ! LOCAL
   INTEGER :: iOp
   
      thisStorage % nDOF = nDOF
      thisStorage % nOps = nOps
      thisStorage % nPeriods = nPeriods
      thisStorage % nTracers = nTracers
      thisStorage % opPeriod = opPeriod
      thisStorage % dt = dt
      thisStorage % CurrentPeriod = -1
   
      ALLOCATE( thisStorage % transportOps(1:nOps) )
      
      DO iOp = 1, nOps
         CALL thisStorage % transportOps(iOp) % Build( nDOF, nDOF, nElems(iOp) )
      ENDDO

      ALLOCATE( thisStorage % tracers(1:nDOF,1:nTracers) )
      ALLOCATE( thisStorage % source(1:nDOF,1:nTracers) )
      ALLOCATE( thisStorage % rFac(1:nDOF,1:nTracers) )
      ALLOCATE( thisStorage % mask(1:nDOF,1:nTracers) )
      ALLOCATE( thisStorage % hSetTracers(1:nDOF,1:nTracers) )

      ! Initiallize values to zero
      thisStorage % tracers     = 0.0_prec
      thisStorage % source      = 0.0_prec
      thisStorage % rFac        = 0.0_prec
      thisStorage % mask        = 1.0_prec
      thisStorage % hSetTracers = 0.0_prec

 END SUBROUTINE Build_TracerStorage
! 
 SUBROUTINE Trash_TracerStorage( thisStorage )
 ! S/R Trash
 !
 !    Deallocates memory for the tracer storage class.
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( TracerStorage ), INTENT(inout) :: thisStorage
   ! LOCAL
   INTEGER :: iOp
   
      DO iOp = 1, thisStorage % nOps
         CALL thisStorage % transportOps(iOp) % Trash( )
      ENDDO
      
      DEALLOCATE( thisStorage % transportOps )
      DEALLOCATE( thisStorage % tracers )
      DEALLOCATE( thisStorage % source )
      DEALLOCATE( thisStorage % rFac )
      DEALLOCATE( thisStorage % mask )
      DEALLOCATE( thisStorage % hSetTracers )

 END SUBROUTINE Trash_TracerStorage
!
!
!==================================================================================================!
!----------------------------------------- Type Specific ------------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE LoadSparseConnectivities_TracerStorage( thisStorage, filebase )
 ! 
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( TracerStorage ), INTENT(inout) :: thisStorage
   CHARACTER(*)                          :: fileBase
   ! Local
   CHARACTER(5) :: periodChar


      WRITE( periodChar, '(I5.5)' ) 1
      PRINT*, '  Loading Sparse Connectivity : '//TRIM(fileBase)//'_advect.'//periodChar
      CALL thisStorage % transportOps(1) % ReadSparseConnectivity( TRIM(fileBase)//'_advect.'//periodChar ) 
      PRINT*, '  Loading Sparse Connectivity : '//TRIM(fileBase)//'_vdiffu.'//periodChar
      CALL thisStorage % transportOps(2) % ReadSparseConnectivity( TRIM(fileBase)//'_vdiffu.'//periodChar ) 
      PRINT*, '  Done!'
 
 END SUBROUTINE LoadSparseConnectivities_TracerStorage
!
 SUBROUTINE CheckForNewOperator_TracerStorage( thisStorage, t, filebase )
 ! 
 ! This routine takes in the simulation time "t" and determines which operator
 ! period needs to be loaded. If the operator period changes, then a new
 ! operator is loaded.
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( TracerStorage ), INTENT(inout) :: thisStorage
   REAL(prec), INTENT(in)                :: t
   CHARACTER(*)                          :: fileBase
   ! Local
   INTEGER      :: nthCycle, newPeriod
   REAL(prec)   :: adjT
   CHARACTER(5) :: periodChar

      IF( thisStorage % nPeriods == 1 )THEN
         newPeriod = 0
      ELSE
         nthCycle = INT( t/(REAL(thisStorage % nPeriods,prec)*thisStorage % opPeriod) )
         adjT     = t - REAL(nthCycle*thisStorage % nPeriods)*thisStorage % opPeriod

         newPeriod = INT( adjT/thisStorage % opPeriod )

         IF( adjT == 0.0_prec )THEN
            newPeriod = 0
         ENDIF
      ENDIF

      IF( newPeriod /= thisStorage % currentPeriod )THEN
         PRINT*, t, adjT, nthCycle, newPeriod
        ! IF( thisStorage % currentPeriod == -1 )THEN
        !    WRITE( periodChar, '(I5.5)' ) newPeriod
        ! ELSE
            WRITE( periodChar, '(I5.5)' ) newPeriod+1
        ! ENDIF

         PRINT*, '  Loading Operator : '//TRIM(fileBase)//'_advect.'//periodChar
         CALL thisStorage % transportOps(1) % ReadSparseConnectivity( TRIM(fileBase)//'_advect.'//periodChar ) 
         CALL thisStorage % transportOps(1) % ReadMatrixData( TRIM(fileBase)//'_advect.'//periodChar ) 
         PRINT*, '  Loading Operator : '//TRIM(fileBase)//'_vdiffu.'//periodChar
         CALL thisStorage % transportOps(2) % ReadSparseConnectivity( TRIM(fileBase)//'_vdiffu.'//periodChar ) 
         CALL thisStorage % transportOps(2) % ReadMatrixData( TRIM(fileBase)//'_vdiffu.'//periodChar ) 

         IF( thisStorage % nPeriods == 1 )THEN
            thisStorage % currentPeriod = 0
         ELSE 
            thisStorage % currentPeriod = newPeriod
         ENDIF

      PRINT*, '  Done!'
      ENDIF
 
 END SUBROUTINE CheckForNewOperator_TracerStorage
!  
 SUBROUTINE MaskField_TracerStorage( thisStorage, i, field )
 ! S/R MaskField
 !
 !   Applies the i-th mask to the "field" array. This routine assumes that mask(j) = 1 implies
 !   that the j-th tracer point is a degree of freedom, ie, it is not hard set.
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( TracerStorage ), INTENT(in) :: thisStorage
   INTEGER, INTENT(in)                :: i
   REAL(prec), INTENT(inout)          :: field(1:thisStorage % nDOF)
   
      field = field*thisStorage % mask(:,i)

 END SUBROUTINE MaskField_TracerStorage
!
!
!
 SUBROUTINE MaskTracers_TracerStorage( thisStorage )
 ! S/R MaskTracers
 !
 !    This subroutine masks the tracers and sets the mask-points to the hard-set values specified
 !    by the attribute "hSetTracers". It is assumed that the hard-set tracer field is set to 
 !    zero where the tracer field is not masked out.
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( TracerStorage ), INTENT(inout) :: thisStorage
   ! LOCAL
   INTEGER    :: i, j
   REAL(prec) :: field(1:thisStorage % nDOF)
   
      DO i = 1, thisStorage % nTracers

         field = thisStorage % tracers(:,i)
         CALL thisStorage % MaskField( i, field )

         DO j = 1, thisStorage % nDOF
            field(j) = field(j) + thisStorage % hSetTracers(j,i)
         ENDDO

         thisStorage % tracers(:,i)  = field
 
      ENDDO

 END SUBROUTINE MaskTracers_TracerStorage
!
 SUBROUTINE CalculateTendency_TracerStorage( thisStorage, tracerfield, t, modelflag, tendency )
 ! CalculateTendency
 !
 !  This subroutine manages the call to the correct model for calculating the tendency of a 
 !  passive tracer field
 !  
 !  *Model flag is passed in as a way to choose from a variety of model options. 
 !   See the code documentation for adding a new tracer model.
 ! 
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( TracerStorage ), INTENT(in) :: thisStorage
   REAL(prec), INTENT(in)             :: tracerfield(1:thisStorage % nDOF, 1:thisStorage % nTracers)
   REAL(prec), INTENT(out)            :: tendency(1:thisStorage % nDOF, 1:thisStorage % nTracers)
   REAL(prec), INTENT(in)             :: t
   INTEGER, INTENT(in)                :: modelflag
   ! Local
   INTEGER :: itracer, i


      IF( modelflag == DyeModel )THEN
         tendency = PassiveDyeModel( thisStorage % nOps, &
                                     thisStorage % nTracers, &
                                     thisStorage % nDOF, &
                                     thisStorage % transportOps, &
                                     thisStorage % source, &
                                     thisStorage % rFac, &
                                     tracerField)
      ELSEIF( modelflag == RadioNuclideModel )THEN

         tendency = ParticulateRadioNuclideModel( thisStorage % nOps, &
                                     thisStorage % nTracers, &
                                     thisStorage % nDOF, &
                                     thisStorage % transportOps, &
                                     thisStorage % source, &
                                     thisStorage % rFac, &
                                     tracerField)

      ELSEIF( modelflag == SettlingModel )THEN

         tendency = ParticulateSettlingModel( thisStorage % nOps, &
                                     thisStorage % nTracers, &
                                     thisStorage % nDOF, &
                                     thisStorage % transportOps, &
                                     thisStorage % source, &
                                     thisStorage % rFac, &
                                     tracerField)

      ELSE
         PRINT*,' Module TracerStorage_Class.f90 : S/R CalculateTendency : Unknown Model Flag', modelflag
         PRINT*,' Stopping! '
         STOP
      ENDIF

      tendency(:,thisStorage % nTracers) = VolumeCorrectionTendency( thisStorage % nDOF, thisStorage % transportOps(1) )
      
      DO itracer = 1, thisStorage % nTracers
         DO i = 1, thisStorage % nDOF
            IF( ABS(thisStorage % mask(i,itracer)) > 1.0_prec )THEN
               PRINT*, 'BUST! Improper Mask Value', thisStorage % mask(i,itracer), i, itracer
            ENDIF
            tendency(i,itracer) = tendency(i,itracer)*thisStorage % mask(i,itracer)
         ENDDO
      ENDDO


 END SUBROUTINE CalculateTendency_TracerStorage
!
 FUNCTION PassiveDyeModel( nOperators, nTracers, nDOF,transportOperators, source, rfac, tracerfield ) RESULT( tendency )
 ! PassiveDyeModel
 !
 !  This function calculates the tendency at time t via a matrix vector multiplication of the 
 !  interpolated transport operator and the tracer field.
 !  
 !  *Model flag is passed in as a way to choose from a variety of model options. In this base code
 ! 
 ! =============================================================================================== !
   IMPLICIT NONE
   INTEGER         :: nOperators, nTracers, nDOF
   TYPE(CRSMatrix) :: transportOperators(1:nOperators)
   REAL(prec)      :: tracerfield(1:nDOF, 1:nTracers)
   REAL(prec)      :: source(1:nDOF, 1:nTracers)
   REAL(prec)      :: rfac(1:nDOF, 1:nTracers)
   REAL(prec)      :: tendency(1:nDOF, 1:nTracers)
   ! LOCAL
   INTEGER         :: itracer, iOp

      ! Calculate the contribution from the transport operator
      DO itracer = 1, nTracers-1 ! Only the passive tracers

         ! Advect the tracer (vertical diffusion is done implicitly)
         tendency(1:nDOF,itracer) = transportOperators(1) % MatVecMul( tracerfield(1:nDOF,iTracer) ) 

         ! Add in the relaxation term
         tendency(1:nDOF,itracer) = tendency(1:nDOF,itracer) + &
                                   (source(1:nDOF,itracer) - tracerfield(1:nDOF,iTracer))*&
                                    rFac(1:nDOF,itracer)
      
      ENDDO


 END FUNCTION PassiveDyeModel
!
 FUNCTION ParticulateSettlingModel( nOperators, nTracers, nDOF,transportOperators, source, rfac,tracerfield ) RESULT( tendency )
 ! ParticulateSettlingModel
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   INTEGER         :: nOperators, nTracers, nDOF
   TYPE(CRSMatrix) :: transportOperators(1:nOperators)
   REAL(prec)      :: tracerfield(1:nDOF, 1:nTracers)
   REAL(prec)      :: source(1:nDOF, 1:nTracers)
   REAL(prec)      :: rfac(1:nDOF, 1:nTracers)
   REAL(prec)      :: tendency(1:nDOF, 1:nTracers)
   ! LOCAL
   INTEGER         :: itracer, iOp, iEl
   REAL(prec)      :: p, k, sWeight

      DO iTracer = 1, ntracers-1
         ! Advect and settle the particulate field
         tendency(:,itracer) = transportOperators(1) % MatVecMul( tracerfield(:,itracer) ) 
         tendency(:,itracer) = tendency(:,itracer) + transportOperators(3) % MatVecMul( tracerfield(:,itracer) ) 

         tendency(1:nDOF,itracer) = tendency(1:nDOF,itracer) + &
                                   (source(1:nDOF,itracer) - tracerfield(1:nDOF,iTracer))*&
                                    rFac(1:nDOF,itracer)
      ENDDO

      
 END FUNCTION ParticulateSettlingModel
!
 FUNCTION ParticulateRadioNuclideModel( nOperators, nTracers, nDOF, transportOperators, source, rfac, tracerfield ) RESULT( tendency )
 ! 
 ! =============================================================================================== !
   IMPLICIT NONE
   INTEGER         :: nOperators, nTracers, nDOF
   TYPE(CRSMatrix) :: transportOperators(1:nOperators)
   REAL(prec)      :: tracerfield(1:nDOF, 1:nTracers)
   REAL(prec)      :: source(1:nDOF, 1:nTracers)
   REAL(prec)      :: rfac(1:nDOF, 1:nTracers)
   REAL(prec)      :: tendency(1:nDOF, 1:nTracers)
   ! LOCAL
   INTEGER         :: itracer, iOp, iEl
   REAL(prec)      :: p, k, sWeight

      DO iTracer = 1, ntracers-1
         ! Perform the advection and diffusion of the particulate field
         tendency(:,itracer) = transportOperators(1) % MatVecMul( tracerfield(:,itracer) ) 
         tendency(:,itracer) = tendency(:,itracer) + transportOperators(2) % MatVecMul( tracerfield(:,itracer) ) 
      ENDDO

      ! Perform the particulate settling -> Assumes the second transport operator is
      ! the settling operator, and the first field is the particulate field
      tendency(:,1) = tendency(:,1) +&
                      transportOperators(3) % MatVecMul( tracerfield(:,1) )

      ! Build the "scavenging operator"
      ! This operator is a weighted version of the settling operator.
      ! Each column of the scavenging operator is multiplied by 
      !  kP/(1+kP) 
      ! where "k" is a scavenging efficiency coefficient, and "P" is the corresponding 
      ! particulate field.
      ! The scavenging coefficient is stored in the "rFac" field
      ! Notice that the 3rd transport operator houses the scavenging operator
      DO iEl = 1, transportOperators(4) % nElems

         p   = tracerfield(transportOperators(4) % col(iEl),1) !obtain the particulate concentration
         k   = rfac(transportOperators(4) % col(iEl), 1) ! and the scavenging coefficient
         sWeight = k*p/( 1.0_prec + k*p ) ! Calculate the equilibrium scavenging flux per unit radionuclide
         transportOperators(4) % A(iEl) = transportOperators(3) % A(iEl)*sWeight

      ENDDO

      ! Add in the scavenging operator for the radionuclide
      tendency(:,2) = tendency(:,2) + transportOperators(4) % MatVecMul( tracerfield(:,2) )
      
 END FUNCTION ParticulateRadioNuclideModel
!
 FUNCTION VolumeCorrectionTendency( nDOF, advectiveOperator ) RESULT( tendency )

   IMPLICIT NONE
   INTEGER           :: nDOF
   TYPE( CRSMatrix ) :: advectiveOperator
   REAL(prec)        :: tendency(1:nDOF)
   ! Local
   REAL(prec)        :: fieldOfOnes(1:nDOF)

      fieldOfOnes = 1.0_prec
      tendency = advectiveOperator % MatVecMul( fieldOfOnes )

 END FUNCTION VolumeCorrectionTendency

END MODULE TracerStorage_Class
