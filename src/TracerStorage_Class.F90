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
#include "FEOTS_Macros.h"

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
      REAL(prec), ALLOCATABLE        :: volume(:,:)
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


 CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE Build_TracerStorage( thisStorage, nDOF, nOps, nElems, nPeriods, nTracers, opPeriod, dt, myRank, nProcs)
#undef __FUNC__
#define __FUNC__ "Build_TracerStorage"
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
      ALLOCATE( thisStorage % volume(1:nDOF,1:thisStorage % nTracers) )
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
#undef __FUNC__
#define __FUNC__ "CheckForNewOperator_TracerStorage"
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
         WRITE( periodChar, '(I5.5)' ) newPeriod+1

        INFO('Loading Operator : '//TRIM(fileBase)//'/transport.'//periodChar//'.h5')
        CALL thisStorage % transportOp % ReadCRSMatrix_HDF5( TRIM(fileBase)//'/transport.'//periodChar//'.h5', myRank, nProc ) 
        INFO('Loading Operator : '//TRIM(fileBase)//'/diffusion.'//periodChar//'.h5')
        CALL thisStorage % diffusionOp % ReadCRSMatrix_HDF5( TRIM(fileBase)//'/diffusion.'//periodChar//'.h5', myRank, nProc ) 
         opSwapped = .TRUE.

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
 SUBROUTINE CalculateTendency_TracerStorage( this, tracerfield, t, modelflag, tendency, volCorrection )
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
   CLASS( TracerStorage ), INTENT(in) :: this
   REAL(prec), INTENT(in)             :: tracerfield(1:this % nDOF, 1:this % nTracers)
   REAL(prec), INTENT(out)            :: tendency(1:this % nDOF, 1:this % nTracers)
   REAL(prec), INTENT(out)            :: volCorrection(1:this % nDOF, 1:this % nTracers)
   REAL(prec), INTENT(in)             :: t
   INTEGER, INTENT(in)                :: modelflag
   ! Local
   INTEGER    :: iTracer, i, row, iel, col
   REAL(prec) :: fieldOfOnes(1:this % nDOF)


      DO iTracer = 1, this % nTracers
         DO row = 1, this % transportOp % nRows

            volCorrection(row,iTracer) = 0.0_prec
            tendency(row,iTracer) = (this % source(row,iTracer) - tracerfield(row,iTracer))*this % rFac(row,iTracer)

            DO iel = this % transportOp % rowBounds(1,row), this % transportOp % rowBounds(2,row)
               col = this % transportOp % col(iel)
               tendency(row,iTracer) = tendency(row,iTracer) + this % transportOp % A(iel)*tracerfield(col,iTracer)
               volCorrection(row,iTracer) = volCorrection(row,iTracer) + this % transportOp % A(iel)
            ENDDO
            tendency(row,iTracer) = tendency(row,iTracer)*this % mask(row,iTracer)
            volCorrection(row,iTracer) = volCorrection(row,iTracer)*this % mask(row,iTracer)
          
         ENDDO
      ENDDO
      
 END SUBROUTINE CalculateTendency_TracerStorage
!
END MODULE TracerStorage_Class
