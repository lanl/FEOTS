! TracerStorage_Class.f90
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
!
! ================================================================================================ !
!

   TYPE TracerStorage 
      INTEGER                        :: nDOF, nPeriods, nOps, nTracers
      INTEGER                        :: currentPeriod
      REAL(prec)                     :: opPeriod, dt
      TYPE( CRSMatrix )              :: transportOp 
      TYPE( CRSMatrix )              :: diffusionOp 
      REAL(prec), ALLOCATABLE        :: tracers(:,:)
      REAL(prec), ALLOCATABLE        :: volume(:)
      REAL(prec), ALLOCATABLE        :: source(:,:), rFac(:,:)
      REAL(prec), ALLOCATABLE        :: mask(:,:)
      INTEGER, ALLOCATABLE           :: tracerIDs(:) ! Tracer ID lower and upper bounds

      CONTAINS
      
      PROCEDURE :: Build => Build_TracerStorage
      PROCEDURE :: Trash => Trash_TracerStorage

      PROCEDURE :: CheckForNewOperator => CheckForNewOperator_TracerStorage
      PROCEDURE :: MaskField         => MaskField_TracerStorage
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
 SUBROUTINE Build_TracerStorage( thisStorage, nDOF, nOps, nElems, nPeriods, nTracers, opPeriod, dt, myRank, nProcs)
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
   INTEGER, INTENT(in)                 :: myRank, nProcs
   ! LOCAL
   INTEGER :: iOp, nT, remainder, i
   
      thisStorage % nDOF = nDOF
      thisStorage % nOps = nOps
      thisStorage % nPeriods = nPeriods
      thisStorage % opPeriod = opPeriod
      thisStorage % dt = dt
      thisStorage % CurrentPeriod = -1

      ! myRank varies from 0 to nProcs-1
      ! We need to share nTracers as evenly as possible across
      ! mpi ranks
      nT = nTracers/nProcs
      remainder = nTracers - nT*nProcs
 
      ! The lower bound for the tracer IDs this rank is responsible for
      IF( myRank == nProcs-1 )THEN
        thisStorage % nTracers = nT+remainder
      ELSE
        thisStorage % nTracers = nT
      ENDIF
      
      ALLOCATE( thisStorage % tracerIDs(1:thisStorage % nTracers) )
      ALLOCATE( thisStorage % tracers(1:nDOF,1:thisStorage % nTracers) )
      ALLOCATE( thisStorage % volume(1:nDOF) )
      ALLOCATE( thisStorage % source(1:nDOF,1:thisStorage % nTracers) )
      ALLOCATE( thisStorage % rFac(1:nDOF,1:thisStorage % nTracers) )
      ALLOCATE( thisStorage % mask(1:nDOF,1:thisStorage % nTracers) )

      DO i = 1, thisStorage % nTracers
        thisStorage % tracerIds(i) =  myRank*nT+i
      ENDDO

      ! Initiallize values to zero
      thisStorage % tracers     = 0.0_prec
      thisStorage % volume      = 0.0_prec
      thisStorage % source      = 0.0_prec
      thisStorage % rFac        = 0.0_prec
      thisStorage % mask        = 1.0_prec

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
   
      DEALLOCATE( thisStorage % tracers )
      DEALLOCATE( thisStorage % volume )
      DEALLOCATE( thisStorage % source )
      DEALLOCATE( thisStorage % rFac )
      DEALLOCATE( thisStorage % mask )

 END SUBROUTINE Trash_TracerStorage
!
!
!==================================================================================================!
!----------------------------------------- Type Specific ------------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE CheckForNewOperator_TracerStorage( thisStorage, t, filebase, OpSwapped, myRank, nProc )
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
   LOGICAL                               :: OpSwapped
   INTEGER                               :: myRank, nProc
   ! Local
   INTEGER      :: nthCycle, newPeriod
   REAL(prec)   :: adjT
   CHARACTER(5) :: periodChar

      opSwapped = .FALSE.
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
         WRITE( periodChar, '(I5.5)' ) newPeriod+1

        PRINT*, '  Loading Operator : '//TRIM(fileBase)//'/transport.'//periodChar//'.h5'
        CALL thisStorage % transportOp % ReadCRSMatrix_HDF5( TRIM(fileBase)//'/transport.'//periodChar//'.h5', myRank, nProc ) 
        PRINT*, '  Loading Operator : '//TRIM(fileBase)//'/diffusion.'//periodChar//'.h5'
        CALL thisStorage % diffusionOp % ReadCRSMatrix_HDF5( TRIM(fileBase)//'/diffusion.'//periodChar//'.h5', myRank, nProc ) 
         opSwapped = .TRUE.

        PRINT*, 'myRank, min/max transport operator : ',myRank, &
                                                        MINVAL(thisStorage % transportOp % A), &
                                                        MAXVAL(thisStorage % transportOp % A)
         IF( thisStorage % nPeriods == 1 )THEN
            thisStorage % currentPeriod = 0
         ELSE 
            thisStorage % currentPeriod = newPeriod
         ENDIF

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
 SUBROUTINE CalculateTendency_TracerStorage( thisStorage, tracerfield, t, modelflag, tendency, volCorrection )
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
   REAL(prec), INTENT(out)            :: volCorrection(1:thisStorage % nDOF)
   REAL(prec), INTENT(in)             :: t
   INTEGER, INTENT(in)                :: modelflag
   ! Local
   INTEGER    :: iTracer, i, row, iel
   REAL(prec) :: fieldOfOnes(1:thisStorage % nDOF)


!      PRINT*, 'Pre-Mask : Tracer, Tracer min/max : ', thisStorage % tracerIDs(1),&
!                                             MINVAL(tracerfield(:,1)),&
!                                             MAXVAL(tracerfield(:,1))
      tendency = PassiveDyeModel( thisStorage % nOps, &
                                  thisStorage % nTracers, &
                                  thisStorage % nDOF, &
                                  thisStorage % transportOp, &
                                  thisStorage % source, &
                                  thisStorage % rFac, &
                                  tracerField)
!      PRINT*, 'Pre-Mask : Tracer, Tendency min/max : ', thisStorage % tracerIDs(1),&
!                                             MINVAL(tendency(:,1)),&
!                                             MAXVAL(tendency(:,1))
!

      DO iTracer = 1, thisStorage % nTracers
         !$OMP DO
         DO i = 1, thisStorage % nDOF
            tendency(i,iTracer) = tendency(i,iTracer)*thisStorage % mask(i,iTracer)
         ENDDO
         !$OMP ENDDO
      ENDDO

      PRINT*, 'Post-Mask : Tracer, Tendency min/max : ', thisStorage % tracerIDs(1),&
                                             MINVAL(tendency(:,1)),&
                                             MAXVAL(tendency(:,1))

      ! Calculate the cell volume update
      !$OMP DO  
      DO row = 1, thisStorage % transportOp % nRows
         volCorrection(row) = 0.0_prec
         DO iel = thisStorage % transportOp % rowBounds(1,row), thisStorage % transportOp % rowBounds(2,row)
            volCorrection(row) = volCorrection(row) + thisStorage % transportOp % A(iel)
         ENDDO
      ENDDO
      !$OMP END DO

      DO iTracer = 1, thisStorage % nTracers
        !$OMP DO  
        DO i = 1, thisStorage % nDOF
         volcorrection(i) = volcorrection(i)*thisStorage % mask(i,iTracer)
        ENDDO
        !$OMP END DO
      ENDDO

 END SUBROUTINE CalculateTendency_TracerStorage
!
 FUNCTION DiffusiveAction( diffusiveOperator, tracerField, nTracers, nDOF ) RESULT( tendency )
   IMPLICIT NONE
   TYPE(CRSMatrix) :: diffusiveOperator
   INTEGER         :: nDOF, nTracers
   REAL(prec)      :: tracerField(1:nDOF,1:nTracers)
   REAL(prec)      :: tendency(1:nDOF,1:nTracers)
   ! Local 
   INTEGER :: i, iTracer, row, iEl, col

     DO iTracer = 1, nTracers ! Only the passive tracers

        !$OMP DO
        DO row = 1, diffusiveOperator % nRows
           tendency(row,iTracer) = 0.0_prec
           DO iel = diffusiveOperator % rowBounds(1,row), diffusiveOperator % rowBounds(2,row)
              col = diffusiveOperator % col(iel)
              tendency(row,iTracer) = tendency(row,iTracer) + diffusiveOperator % A(iel)*tracerfield(col,iTracer)
           ENDDO
        ENDDO
        !$OMP ENDDO

     ENDDO

 END FUNCTION DiffusiveAction
!
 FUNCTION PassiveDyeModel( nOperators, nTracers, nDOF,transportOperator, source, rfac, tracerfield ) RESULT( tendency )
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
   TYPE(CRSMatrix) :: transportOperator
   REAL(prec)      :: tracerfield(1:nDOF, 1:nTracers)
   REAL(prec)      :: source(1:nDOF, 1:nTracers)
   REAL(prec)      :: rfac(1:nDOF, 1:nTracers)
   REAL(prec)      :: tendency(1:nDOF, 1:nTracers)
   ! LOCAL
   INTEGER         :: iTracer, i, row, col, iel

      ! Calculate the contribution from the transport operator
      DO iTracer = 1, nTracers ! Only the passive tracers

         ! Advect the tracer (vertical diffusion is done implicitly)
         !$OMP DO
         DO row = 1, transportOperator % nRows
            tendency(row,iTracer) = 0.0_prec
            DO iel = transportOperator % rowBounds(1,row), transportOperator % rowBounds(2,row)
               col = transportOperator % col(iel)
               tendency(row,iTracer) = tendency(row,iTracer) + transportOperator % A(iel)*tracerfield(col,iTracer)
            ENDDO
         ENDDO
         !$OMP ENDDO

         ! Add in the relaxation term
         !$OMP DO
         DO i = 1, nDOF
            tendency(i,iTracer) = tendency(i,iTracer) + &
                                      (source(i,iTracer) - tracerfield(i,iTracer))*&
                                       rFac(i,iTracer)
         ENDDO
         !$OMP ENDDO

      ENDDO
     
 END FUNCTION PassiveDyeModel
!
! FUNCTION ParticulateSettlingModel( nOperators, nTracers, nDOF,transportOperators, source, rfac,tracerfield ) RESULT( tendency )
! ! ParticulateSettlingModel
! !
! ! =============================================================================================== !
!   IMPLICIT NONE
!   INTEGER         :: nOperators, nTracers, nDOF
!   TYPE(CRSMatrix) :: transportOperators(1:nOperators)
!   REAL(prec)      :: tracerfield(1:nDOF, 1:nTracers)
!   REAL(prec)      :: source(1:nDOF, 1:nTracers)
!   REAL(prec)      :: rfac(1:nDOF, 1:nTracers)
!   REAL(prec)      :: tendency(1:nDOF, 1:nTracers)
!   ! LOCAL
!   INTEGER         :: iTracer, iOp, iEl
!   REAL(prec)      :: p, k, sWeight
!
!      DO iTracer = 1, ntracers
!         ! Advect and settle the particulate field
!
!         !$OMP PARALLEL
!         tendency(:,iTracer) = transportOperators(1) % MatVecMul( tracerfield(:,iTracer) ) 
!         !$OMP ENDPARALLEL
!
!         !$OMP PARALLEL
!         tendency(:,iTracer) = tendency(:,iTracer) + transportOperators(3) % MatVecMul( tracerfield(:,iTracer) ) 
!         !$OMP ENDPARALLEL
!
!         tendency(1:nDOF,iTracer) = tendency(1:nDOF,iTracer) + &
!                                   (source(1:nDOF,iTracer) - tracerfield(1:nDOF,iTracer))*&
!                                    rFac(1:nDOF,iTracer)
!      ENDDO
!
!      
! END FUNCTION ParticulateSettlingModel
!!
! FUNCTION ParticulateRadioNuclideModel( nOperators, nTracers, nDOF, transportOperators, source, rfac, tracerfield ) RESULT( tendency )
! ! 
! ! =============================================================================================== !
!   IMPLICIT NONE
!   INTEGER         :: nOperators, nTracers, nDOF
!   TYPE(CRSMatrix) :: transportOperators(1:nOperators)
!   REAL(prec)      :: tracerfield(1:nDOF, 1:nTracers)
!   REAL(prec)      :: source(1:nDOF, 1:nTracers)
!   REAL(prec)      :: rfac(1:nDOF, 1:nTracers)
!   REAL(prec)      :: tendency(1:nDOF, 1:nTracers)
!   ! LOCAL
!   INTEGER         :: iTracer, iOp, iEl
!   REAL(prec)      :: p, k, sWeight
!
!      DO iTracer = 1, ntracers
!         ! Perform the advection and diffusion of the particulate field
!         tendency(:,iTracer) = transportOperators(1) % MatVecMul( tracerfield(:,iTracer) ) 
!      !   tendency(:,iTracer) = tendency(:,iTracer) + transportOperators(2) % MatVecMul( tracerfield(:,iTracer) ) 
!      ENDDO
!
!      ! Perform the particulate settling -> Assumes the second transport operator is
!      ! the settling operator, and the first field is the particulate field
!      tendency(:,1) = tendency(:,1) +&
!                      transportOperators(3) % MatVecMul( tracerfield(:,1) )
!
!      ! Build the "scavenging operator"
!      ! This operator is a weighted version of the settling operator.
!      ! Each column of the scavenging operator is multiplied by 
!      !  kP/(1+kP) 
!      ! where "k" is a scavenging efficiency coefficient, and "P" is the corresponding 
!      ! particulate field.
!      ! The scavenging coefficient is stored in the "rFac" field
!      ! Notice that the 3rd transport operator houses the scavenging operator
!      DO iEl = 1, transportOperators(4) % nElems
!
!         p   = tracerfield(transportOperators(4) % col(iEl),1) !obtain the particulate concentration
!         k   = rfac(transportOperators(4) % col(iEl), 1) ! and the scavenging coefficient
!         sWeight = k*p/( 1.0_prec + k*p ) ! Calculate the equilibrium scavenging flux per unit radionuclide
!         transportOperators(4) % A(iEl) = transportOperators(3) % A(iEl)*sWeight
!
!      ENDDO
!
!      ! Add in the scavenging operator for the radionuclide
!      tendency(:,2) = tendency(:,2) + transportOperators(4) % MatVecMul( tracerfield(:,2) )
!      ! Add in the uniform source for the radionuclide
!      tendency(:,2) = tendency(:,2) + source(1:nDOF,2)
!      
! END FUNCTION ParticulateRadioNuclideModel
!
END MODULE TracerStorage_Class
