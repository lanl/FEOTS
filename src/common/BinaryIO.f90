! BinaryIO.f90
! 
! Copyright 2016 Joseph Schoonover, Los Alamos National Laboratory (jschoonover@lanl.gov) 
! All rights reserved. 
! 
! BinaryIO.f90 is part of the Fast Equilibration of Ocean Tracers Software (FEOTS). 
! 
! Redistribution and use in source and binary forms, with or without modification, are permitted 
! provided that the following conditions are met: 
! 
! 1. Redistributions of source code must retain the above copyright notice, this list of conditions 
! and the following disclaimer. 
! 
! 2. Redistributions in binary form must reproduce the above copyright notice, this list of 
! conditions and the following disclaimer in the documentation and/or other materials provided with 
! the distribution. 
! 
! 3. Neither the name of the copyright holder nor the names of its contributors may be used to 
! endorse or promote products derived from this software without specific prior written permission. 
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR  
! IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND  
! FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR 
! CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
! WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY 
! WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !


MODULE BinaryIO
!
! BinaryIO.f90 ( v 1.1 )
!
! This module was written to provide subroutines which READ binary files that are
! written in big-endian byte ordering. The binary files can be READ into 1-,2-, or
! 3-Dimensional arrays.
!
! >>>>>>>>>>>>> List of Functions and Subroutines <<<<<<<<<<<<<<
!
!
! o  ReadBIN(1D/2D/3D) ------------------------------------------------------------------------------- 
!
!      The ReadBin routines WRITE a binary FILE with big-endian byte ordering. The length of the FILE
!      is determined by the precision and the total number of elements in an array. The only required
!      input is the filename and the variable with the appropriate size. An flag is returned to help
!      the outside calling program determine what to DO next, in the event of an unsuccessful READ.
!
!       ** The ReadBIN(1D/2D/3D) routines are overloaded onto the one procedure "ReadBIN" so that the
!          user should only have to call ReadBIN( filename, variable, ioFlag ) without the extension of 
!          "1D", "2D", or "3D". 
!
! o  WriteBIN(1D/2D/3D) ------------------------------------------------------------------------------- 
!
!      The WriteBin routines WRITE a binary FILE with big-endian byte ordering. The length of the FILE
!      is determined by the precision and the total number of elements in an array. The only required
!      input is the filename and the variable. 
!
!       ** The WriteBIN(1D/2D/3D) routines are overloaded onto the one procedure "WriteBIN" so that the
!          user should only have to call WriteBIN( filename, variable ) without the extension of 
!          "1D", "2D", or "3D". 
! 
! ========================================================================================================= !

 USE ModelPrecision
 USE CommonRoutines

 IMPLICIT NONE


    ! Overload the subroutine name "ReadBIN" to call the 1,2,or 3-D implementations based
    ! on the array dimensions passed. This way, the user only has to call ReadBIN( theFile, var, fileLength )
    ! to READ in a FILE to a REAL array.
    INTERFACE ReadBIN
        MODULE PROCEDURE ReadBIN4D, ReadBIN3D, ReadBIN2D, ReadBIN1D, ReadBIN2D_INT
    END INTERFACE

    INTERFACE WriteBIN
        MODULE PROCEDURE WriteBIN4D, WriteBIN3D, WriteBIN2D, WriteBIN1D
    END INTERFACE
 

 CONTAINS

! ====================================================================================== !
! ------------------------------- MITgcm Standard IO
! ====================================================================================== !
!
!
!
 SUBROUTINE ReadBIN4D(theFile, var, ioFlag )
 ! S/R ReadBIN4D
 !
 !=============================================================!
 ! Declarations
   IMPLICIT NONE
   CHARACTER(*), INTENT(in) :: theFile
   REAL(prec), INTENT(out)  :: var(1:,1:,1:,1:)
   INTEGER, INTENT(out)     :: ioFlag
   ! Local
   INTEGER :: fUnit, fileLength

     fileLength = SIZE( var )

     OPEN(UNIT=NewUnit(fUnit),FILE=trim(theFile),   &
                              form='unformatted',   &
                              access='direct',      &
                              CONVERT='big_endian', &
                              recl=prec*fileLength )

     READ(UNIT=fUnit,rec=1, iostat=ioFlag) var

     CLOSE(UNIT=fUnit)


 END SUBROUTINE ReadBIN4D
!
!
!
 SUBROUTINE ReadBIN3D(theFile, var, ioFlag )
 ! S/R ReadBIN3D
 !
 !=============================================================!
 ! Declarations
   IMPLICIT NONE
   CHARACTER(*), INTENT(in) :: theFile
   REAL(prec), INTENT(out)  :: var(1:,1:,1:)
   INTEGER, INTENT(out)     :: ioFlag
   ! Local
   INTEGER :: fUnit, fileLength
   

     fileLength = SIZE( var )

     OPEN(UNIT=NewUnit(fUnit),FILE=trim(theFile),   &
                              form='unformatted',   &
                              access='direct',      &
                              CONVERT='big_endian', &
                              recl=prec*fileLength )

     READ(UNIT=fUnit,rec=1, iostat=ioFlag) var

     CLOSE(UNIT=fUnit)


 END SUBROUTINE ReadBIN3D
!
!
!
 SUBROUTINE ReadBIN2D(theFile, var, ioFlag )
 ! S/R ReadBIN2D
 !
 !=============================================================!
 ! Declarations
   IMPLICIT NONE
   CHARACTER(*), INTENT(in) :: theFile
   REAL(prec), INTENT(out)  :: var(1:,1:)
   INTEGER, INTENT(out)     :: ioFlag
   ! Local
   INTEGER :: fUnit, fileLength

     fileLength = SIZE( var )

     OPEN(UNIT=NewUnit(fUnit),FILE=trim(theFile),   &
                              form='unformatted',   &
                              access='direct',      &
                              CONVERT='big_endian', &
                              recl=prec*fileLength )

     READ(UNIT=fUnit,rec=1, iostat=ioFlag) var

     CLOSE(UNIT=fUnit)


 END SUBROUTINE ReadBIN2D
!
!
!
 SUBROUTINE ReadBIN1D(theFile, var, ioFlag )
 ! S/R ReadBIN1D
 !
 !=============================================================!
 ! Declarations
   IMPLICIT NONE
   CHARACTER(*), INTENT(in) :: theFile
   REAL(prec), INTENT(out)  :: var(1:)
   INTEGER, INTENT(out)     :: ioFlag
   ! Local
   INTEGER :: fUnit, fileLength

     fileLength = SIZE( var )

     OPEN(UNIT=NewUnit(fUnit),FILE=trim(theFile),   &
                              form='unformatted',   &
                              access='direct',      &
                              CONVERT='big_endian', &
                              recl=prec*fileLength )

     READ(UNIT=fUnit,rec=1,iostat=ioFlag) var

     CLOSE(UNIT=fUnit)


 END SUBROUTINE ReadBIN1D
!
!
!
 SUBROUTINE ReadBIN2D_INT(theFile, var, ioFlag )
 ! S/R ReadBIN2D
 !
 !=============================================================!
 ! Declarations
   IMPLICIT NONE
   CHARACTER(*), INTENT(in) :: theFile
   INTEGER, INTENT(out)     :: var(1:,1:)
   INTEGER, INTENT(out)     :: ioFlag
   ! Local
   INTEGER :: fUnit, fileLength

     fileLength = SIZEOF( var )

     OPEN(UNIT=NewUnit(fUnit),FILE=trim(theFile),   &
                              form='unformatted',   &
                              access='direct',      &
                              CONVERT='big_endian', &
                              recl=fileLength )

     READ(UNIT=fUnit,rec=1, iostat=ioFlag) var

     CLOSE(UNIT=fUnit)


 END SUBROUTINE ReadBIN2D_INT
!
!
!

 SUBROUTINE WriteBIN4D(theFile, var )
 ! S/R WriteBIN4D
 !
 !=============================================================!
 ! Declarations
   IMPLICIT NONE
   CHARACTER(*), INTENT(in) :: theFile
   REAL(prec), INTENT(in)   :: var(1:,1:,1:,1:)

   ! Local
   INTEGER :: fUnit, fileLength

     fileLength = SIZE( var )

     OPEN(UNIT=NewUnit(fUnit),FILE=trim(theFile),   &
                              form='unformatted',   &
                              access='direct',      &
                              CONVERT='big_endian', &
                              status='replace',     &
                              recl=prec*fileLength )

     WRITE(UNIT=fUnit,rec=1) var

     CLOSE(UNIT=fUnit)


 END SUBROUTINE WriteBIN4D
!
!
!
 SUBROUTINE WriteBIN3D(theFile, var )
 ! S/R WriteBIN3D
 !
 !=============================================================!
 ! Declarations
   IMPLICIT NONE
   CHARACTER(*), INTENT(in) :: theFile
   REAL(prec), INTENT(in)   :: var(1:,1:,1:)

   ! Local
   INTEGER :: fUnit, fileLength

     fileLength = SIZE( var )

     OPEN(UNIT=NewUnit(fUnit),FILE=trim(theFile),   &
                              form='unformatted',   &
                              access='direct',      &
                              CONVERT='big_endian', &
                              status='replace',     &
                              recl=prec*fileLength )

     WRITE(UNIT=fUnit,rec=1) var

     CLOSE(UNIT=fUnit)


 END SUBROUTINE WriteBIN3D
!
!
!
 SUBROUTINE WriteBIN2D(theFile, var )
 ! S/R ReadBIN2D
 !
 !=============================================================!
 ! Declarations
   IMPLICIT NONE
   CHARACTER(*), INTENT(in) :: theFile
   REAL(prec), INTENT(in)   :: var(1:,1:)

   ! Local
   INTEGER :: fUnit, fileLength

     fileLength = SIZE( var )

     OPEN(UNIT=NewUnit(fUnit),FILE=trim(theFile),   &
                              form='unformatted',   &
                              access='direct',      &
                              CONVERT='big_endian', &
                              status='replace',     &
                              recl=prec*fileLength )

     WRITE(UNIT=fUnit,rec=1) var

     CLOSE(UNIT=fUnit)


 END SUBROUTINE WriteBIN2D
!
!
!
 SUBROUTINE WriteBIN1D(theFile, var )
 ! S/R ReadBIN1D
 !
 !=============================================================!
 ! Declarations
   IMPLICIT NONE
   CHARACTER(*), INTENT(in) :: theFile
   REAL(prec), INTENT(in)   :: var(1:)

   ! Local
   INTEGER :: fUnit, fileLength

     fileLength = SIZE( var )

     OPEN(UNIT=NewUnit(fUnit),FILE=trim(theFile),   &
                              form='unformatted',   &
                              access='direct',      &
                              CONVERT='big_endian', &
                              status='replace',     &
                              recl=prec*fileLength )

     WRITE(UNIT=fUnit,rec=1) var

     CLOSE(UNIT=fUnit)


 END SUBROUTINE WriteBIN1D
!
!
!
END MODULE BinaryIO
