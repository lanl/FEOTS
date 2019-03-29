! FEOTSDriver.f90
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

PROGRAM FEOTSDriver


! src/POP/
USE POP_FEOTS_Class
USE POP_Params_Class

#ifdef HAVE_OPENMP
USE OMP_LIB
#endif

IMPLICIT NONE

   TYPE( POP_FEOTS ) :: feots
   CHARACTER(10) :: ncFileTag
   CHARACTER(5)  :: fileIDChar
   CHARACTER(200):: thisIRFFile
   INTEGER       :: funit, recordID, fileID, i, nIODumps
   INTEGER       :: mpiErr, myRank, nProcs, iter
   REAL(prec)    :: tn
   REAL(prec)    :: t1, t2
   
#ifdef HAVE_MPI
      CALL MPI_INIT( mpiErr )
      CALL MPI_COMM_RANK( MPI_COMM_WORLD, myRank, mpiErr )
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD, nProcs, mpiErr )
#else
      myRank = 0
      nProcs = 1
#endif

      CALL feots % Build( myRank, nProcs )

      recordID = 1
      IF( feots % params % runMode == FORWARD )THEN 

         IF( myRank == 0 )THEN 
            fileID   = feots % params % iterInit
            ! /////////////////////////// Load in the initial conditions //////////////////// !
            WRITE( ncfileTag, '(I10.10)' ) fileID
            ! For now, the record ID is 1. In the future, this will need to be
            ! calculated as a function of the initial iterate, the dump frequency,
            ! and the number of records per netcdf file
            
           ! Tracer.init.nc is read for the mask and source terms
            CALL feots % nativeSol % InitializeForNetCDFRead( feots % params % TracerModel,'Tracer.init.nc', .TRUE. )
            CALL feots % nativeSol % ReadSourceEtcNetCDF( feots % mesh )
            CALL feots % nativeSol % ReadNetCDFRecord( feots % mesh, recordID )
            CALL feots % nativeSol % FinalizeNetCDF( )


            ! ************
            ! This section of code loads in the temperature, salinity, potential
            ! density, and ssh fields if the water mass tagging is turned on.
            ! In future implementations with the volume correction, this will need
            ! to be turned on the volume corrections are enabled. 
            !
            IF( feots % params % WaterMassTagging ) THEN

               WRITE( fileIDChar, '(I5.5)' ) feots % params % IRFStart
               IF( feots % params % Regional )THEN
                  CALL feots % nativeSol % LoadOceanState( feots % mesh, &
                                                          TRIM(feots % params % regionalOperatorDirectory)//'Ocean.'//fileIDChar//'.nc')
               ELSE
                  OPEN( UNIT=NewUnit(fUnit),&
                        FILE=TRIM(feots % params % IRFListFile), &
                        FORM='FORMATTED',&
                        ACCESS='SEQUENTIAL',&
                        ACTION='READ',&
                        STATUS='OLD' )
            
                  DO fileID = 1, feots % params % nIRFFiles
            
                     READ( fUnit, '(A200)' ) thisIRFFile
            
                     IF( fileID == feots % params % IRFStart )THEN
                        CALL feots % nativeSol % LoadOceanState( feots % mesh,TRIM(thisIRFFile) )
                     ENDIF
                  ENDDO
   
                  CLOSE(fUnit)
               ENDIF
            ENDIF

         ENDIF ! myRank == 0
         !***********

         IF( feots % params % iterInit == 0 )THEN
            tn = 0.0_prec
         ELSE

            IF( myRank == 0 )THEN 
               ! This pickup file is read for the correct "initial condition"
               CALL feots % nativeSol % InitializeForNetCDFRead( feots % params % TracerModel,'Tracer.'//ncFileTag//'.nc', .FALSE. )
               CALL feots % nativeSol % ReadNetCDFRecord( feots % mesh, recordID )
               CALL feots % nativeSol % FinalizeNetCDF( )
            ENDIF
               tn = REAL(feots % params % iterInit,prec)*feots % params % dt

         ENDIF
   
         ! /////////////////////////////////////////////////////////////////////////////// !
         
         ! Transfer the data from the native storage to the FEOTS storage
         IF( myRank == 0 )THEN
            CALL feots % MapAllToDOF( )
         ENDIF
#ifdef HAVE_MPI
         CALL MPI_BARRIER( MPI_COMM_WORLD, mpiErr )
         CALL feots % ScatterSolution( myRank, nProcs )
         CALL feots % ScatterSource( myRank, nProcs )
         CALL feots % ScatterMask( myRank, nProcs )
         IF( myRank /= 0 )THEN
            CALL feots % MapTracerToDOF( )
         ENDIF
#endif 
         ! //// Forward Mode //// !
         PRINT*, '  Starting ForwardStep'
   
         DO iter = feots % params % iterInit, feots % params % iterInit + feots % params % nTimeSteps -1, feots % params % nStepsPerDump

            !$OMP PARALLEL
            IF( myRank /= 0 .OR. nProcs == 0)THEN
               CALL feots % ForwardStep( tn, feots % params % nStepsPerDump, myRank, nProcs )
            ENDIF
            !$OMP END PARALLEL 

            IF(myRank /= 0  .OR. nProcs == 0)THEN
               CALL feots % MapTracerFromDOF( )
            ENDIF

#ifdef HAVE_MPI
            CALL feots % GatherSolution( myRank, nProcs )
#endif
            IF( myRank == 0 )THEN

               WRITE( ncfileTag, '(I10.10)' ) iter + feots % params % nStepsPerDump
               CALL feots % nativeSol % InitializeForNetCDFWrite( feots % params % TracerModel, &
                                                                  feots % mesh, &
                                                                 'Tracer.'//ncFileTag//'.nc', &
                                                                  .FALSE. )
               CALL feots % nativeSol % WriteNetCDFRecord( feots % mesh, recordID )
               CALL feots % nativeSol % FinalizeNetCDF( )
            ENDIF

         ENDDO
#ifdef HAVE_MPI
         CALL MPI_BARRIER( MPI_COMM_WORLD, mpiErr )
#endif

      ELSEIF( feots % params % runMode == EQUILIBRIUM )THEN

         IF( myRank == 0 )THEN
           ! Tracer.init.nc is read for the mask and source terms
            CALL feots % nativeSol % InitializeForNetCDFRead( feots % params % TracerModel,'Tracer.init.nc', .TRUE. )
            CALL feots % nativeSol % ReadSourceEtcNetCDF( feots % mesh )
            CALL feots % nativeSol % ReadNetCDFRecord( feots % mesh, recordID )
            CALL feots % nativeSol % FinalizeNetCDF( )
   
  
            IF( feots % params % isPickupRun )THEN  
               CALL feots % nativeSol % InitializeForNetCDFRead( feots % params % TracerModel,'Tracer.pickup.nc', .FALSE. )
               CALL feots % nativeSol % ReadNetCDFRecord( feots % mesh, recordID )
               CALL feots % nativeSol % FinalizeNetCDF( )
               tn = REAL(feots % params % iterInit,prec)*feots % params % dt
   
            ENDIF
   
            ! /////////////////////////////////////////////////////////////////////////////// !
            
            ! Transfer the data from the native storage to the FEOTS storage
            CALL feots % MapAllToDOF( )
         ENDIF
#ifdef HAVE_MPI
         CALL MPI_BARRIER( MPI_COMM_WORLD, mpiErr )
         CALL feots % ScatterSolution( myRank, nProcs )
         CALL feots % ScatterSource( myRank, nProcs )
         CALL feots % ScatterMask( myRank, nProcs )
#endif 

#ifdef HAVE_OPENMP   
            t1 = omp_get_wtime( )
#else
            CALL CPU_TIME( t1 )
#endif
         CALL feots % JFNK( myRank )

          
#ifdef HAVE_OPENMP
            t2 = omp_get_wtime( )
#else
            CALL CPU_TIME( t2 )
#endif
            PRINT*, 'JFNK wall time :', t2-t1

      ENDIF

      CALL feots % Trash( )
#ifdef HAVE_MPI
      CALL MPI_FINALIZE( mpiErr )
#endif


END PROGRAM FEOTSDriver
