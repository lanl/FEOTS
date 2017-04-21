! ConstantsDictionary.f90
! 
! Copyright 2016 Joseph Schoonover, Los Alamos National Laboratory (jschoonover@lanl.gov) 
! All rights reserved. 
! 
! //////////////////////////////////////////////////////////////////////////////////////////////// !
 
 

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


 END MODULE ConstantsDictionary
