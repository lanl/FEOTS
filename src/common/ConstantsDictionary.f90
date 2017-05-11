! ConstantsDictionary.f90
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
 
 

 MODULE ConstantsDictionary

  USE ModelPrecision

  
  !*************************************************************!
  ! --------------------- File I/O Behavior --------------------!
  ! ************************************************************!
  INTEGER, PARAMETER :: nInt4PerChunk=25000000 ! 100 MB file IO resolution
  INTEGER, PARAMETER :: nInt8PerChunk=12500000 ! 100 MB file IO resolution
  INTEGER, PARAMETER :: nReal4PerChunk=25000000 ! 100 MB file IO resolution
  INTEGER, PARAMETER :: nReal8PerChunk=12500000 ! 100 MB file IO resolution


  !*************************************************************!
  ! ------------------ MATHEMATICAL CONSTANTS ------------------!
  ! ************************************************************!
  !                                                             !
  ! ------------------------------------------------------------!
   REAL(prec), PARAMETER :: pi   = 4.0_prec*atan(1.0_prec)
   REAL(prec), PARAMETER :: ZERO = 0.0_prec
   REAL(prec), PARAMETER :: ONE  = 1.0_prec
   REAL(prec), PARAMETER :: TWO  = 2.0_prec
   REAL(prec), PARAMETER :: HALF = 0.5_prec
   REAL(prec), PARAMETER :: TOL  = epsilon(1.0_prec)

   INTEGER, PARAMETER :: kItMax = 50 ! Max iterations for Newton's method.
                                      ! USEd in : MODULE LEGENDRE.f90, S/R LEG_GAUSSQUAD
   REAL(prec), PARAMETER :: fillValue = -999.999_prec
  !*************************************************************!
  ! ----------------- ROOT FINDER CONSTANTS --------------------!
  ! ************************************************************!
  !                                                             !
  ! ------------------------------------------------------------!
   REAL(prec), PARAMETER :: newtonTolerance = 10.0**(-10)
   INTEGER, PARAMETER    :: newtonMax       = 100
  
  !*************************************************************!
  ! ----------------- TIME STEPPING CONSTANTS ------------------!
  ! ************************************************************!
  !                                                             !
  ! ------------------------------------------------------------!
  ! Runge-Kutta 3rd Order, low storage constants
   REAL(prec), PARAMETER :: rk3_a(1:3) = (/ 0.0_prec, -5.0_prec/9.0_prec, -153.0_prec/128.0_prec /)
   REAL(prec), PARAMETER :: rk3_b(1:3) = (/ 0.0_prec, 1.0_prec/3.0_prec, 3.0_prec/4.0_prec /)
   REAL(prec), PARAMETER :: rk3_g(1:3) = (/ 1.0_prec/3.0_prec, 15.0_prec/16.0_prec, 8.0_prec/15.0_prec /)

  !*************************************************************!
  ! ------------------- PHYSICAL CONSTANTS ---------------------!
  ! ************************************************************!
  !                                                             !
  ! ------------------------------------------------------------!
  ! Time conversion factors
   REAL(prec), PARAMETER   :: secondsToMinutes = 1.0_prec/60.0_prec                   ! conversion for seconds to minutes
   REAL(prec), PARAMETER   :: minutesToHours   = 1.0_prec/60.0_prec                   ! conversion for minutes to hours
   REAL(prec), PARAMETER   :: hoursToDays      = 1.0_prec/24.0_prec                   ! conversion for hours to days
   REAL(prec), PARAMETER   :: daysToMonths     = 12.0_prec/365.25_prec                ! conversion for days to months
   REAL(prec), PARAMETER   :: monthsToYears    = 1.0_prec/12.0_prec                   ! conversion for months to years
   REAL(prec), PARAMETER   :: daysToSeconds    = 86400.0_prec

  !*************************************************************!
  ! ------------------- Package defaults ---------------------!
  ! ************************************************************!
  !                                                             !
  ! ------------------------------------------------------------!
   INTEGER, PARAMETER :: defaultNameLength = 40

  
  ! Misc. INTEGER and CHARACTER flag definitions
   INTEGER, PARAMETER      :: NONE = 0
   CHARACTER(1), PARAMETER :: nada = ' ' 


  !*************************************************************!
  ! ------------------------- Model Flags ----------------------!
  ! ************************************************************!
   ! Tracer Models !
   INTEGER, PARAMETER :: DyeModel             = 500
   INTEGER, PARAMETER :: RadionuclideModel    = 501
   INTEGER, PARAMETER :: SettlingModel        = 503
   INTEGER, PARAMETER :: ImpulseField         = 502
   INTEGER, PARAMETER :: ImpulseResponseField = 502
   ! Tracer integration modes
   INTEGER, PARAMETER :: forward     = 600
   INTEGER, PARAMETER :: equilibrium = 601

   ! Stencil Flags
   INTEGER, PARAMETER :: LaxWendroff = 700
   INTEGER, PARAMETER :: LaxWendroffOverlap = 701
   
   INTEGER, PARAMETER :: PeriodicTripole = 800
   INTEGER, PARAMETER :: Dipole = 801

   INTEGER, PARAMETER :: Normal  = 850
   INTEGER, PARAMETER :: Overlap = 851
   INTEGER, PARAMETER :: Lateral = 852
   INTEGER, PARAMETER :: LateralPlusCorners = 853
   INTEGER, PARAMETER :: Temperature = 1000
   INTEGER, PARAMETER :: Salinity = 1000
   INTEGER, PARAMETER :: Density = 1000


 END MODULE ConstantsDictionary
