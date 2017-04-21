MODULE POP_GridTypeMappings

! src/common/
USE ConstantsDictionary


IMPLICIT NONE

CONTAINS



 
 SUBROUTINE GetTrueIJ( meshtype, i, j, nX, nY, true_i, true_j )
 
   IMPLICIT NONE
   INTEGER, INTENT(in)    :: meshtype, i, j, nX, nY
   INTEGER, INTENT(inout) :: true_i, true_j

      IF( meshtype == PeriodicTripole )THEN
 
         CALL GetTrueIJ_PeriodicTripole( i, j, nX, nY, true_i, true_j )

      ELSE

         PRINT*, 'Module POP_GridTypeMappings : S/R GetTrueIJ : Unknown Grid Type'
         STOP

      ENDIF

 END SUBROUTINE GetTrueIJ
!
 SUBROUTINE GetTrueIJ_PeriodicTripole( i, j, nX, nY, true_i, true_j )
   IMPLICIT NONE
   INTEGER, INTENT(in)  :: i, j, nX, nY
   INTEGER, INTENT(out) :: true_i, true_j
   ! Local
   INTEGER :: dj

      IF( j > nY )THEN
   
         dj = j - nY
         true_j = nY - (dj-1)
         IF( i == 0 .OR. i == nX+1 )THEN
            true_i = i
         ELSE 
            true_i = nX-(i-1)
         ENDIF
   
      ELSE
   
         true_j = j
          
         ! Taking care of the x-periodicity in the map
         IF( i < 1 )THEN
            true_i = nX + i
         ELSEIF( i > nX )THEN
            true_i = i-nX
         ELSE
            true_i = i
         ENDIF
   
      ENDIF
   
 END SUBROUTINE GetTrueIJ_PeriodicTripole

END MODULE POP_GridTypeMappings
