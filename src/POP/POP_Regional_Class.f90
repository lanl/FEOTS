
MODULE POP_Regional_Class


! src/common/
USE ModelPrecision
! src/POP/
USE POP_Mesh_Class
USE POP_Stencil_Class
USE POP_GridTypeMappings

! This module defines a class for handling masks of a global mesh so that
! regional meshes can be constructed, boundary conditions can be applied, and
! subsets of transport operators extracted


IMPLICIT NONE


! Bounding latitudes should always be specified in degrees North
! Bounding longitudes should always be specified in degrees East

! The "ijk" in ijkInRegion refer to the (i,j,k) triplets from the original
! global mesh. The length of the second dimension of the ijkInRegion array is
! the number of degrees of freedom within the local region.

! The "dof" in dofInRegion refer to the single degree of freedom index from the
! original global mesh. The length of the dofInRegion array is the number of
! degrees of freedom within the local region.
!

! The "boundaryCells" array references local degrees of freedom within the
! regional mesh. Its length is determined by the number of cells that make up
! the border of the regional mesh. 
! To find the physical position of the m-th boundary cell, one can obtain the
! (i,j,k) triplet of the global mesh location by using
!
! (/ i, j, k /) = ijkInRegion(1:3,boundaryCells(m) )
!
!

   TYPE POP_Regional
      REAL(prec) :: south, north, east, west
      LOGICAL    :: crossesPrimeMeridian
      INTEGER    :: nCells, nBCells, nDOF
      INTEGER, ALLOCATABLE :: ijkInRegion(:,:)
      INTEGER, ALLOCATABLE :: dofInRegion(:)  
      INTEGER, ALLOCATABLE :: inverseDOFMap(:)  
      INTEGER, ALLOCATABLE :: boundaryCells(:)
      INTEGER, ALLOCATABLE :: dofToLocalIJK(:,:)

      CONTAINS

      PROCEDURE :: Build => Build_POP_Regional
      PROCEDURE :: Trash => Trash_POP_Regional

      PROCEDURE :: FindThoseInRegion => FindThoseInRegion_POP_Regional 
      PROCEDURE :: FindBoundaryCells => FindBoundaryCells_POP_Regional

      PROCEDURE :: GenerateRegionalMesh => GenerateRegionalMesh_POP_Regional

      PROCEDURE :: WritePickup => WritePickup_POP_Regional
      PROCEDURE :: ReadPickup  => ReadPickup_POP_Regional

   END TYPE POP_Regional


CONTAINS


 SUBROUTINE Build_POP_Regional( myRegion, mesh, mystencil, meshType, south, north, east, west )

   IMPLICIT NONE
   CLASS( POP_Regional ), INTENT(out) :: myRegion
   TYPE( POP_Mesh ), INTENT(inout)    :: mesh
   TYPE( Stencil ), INTENT(in)        :: mystencil
   INTEGER, INTENT(in)                :: meshType
   REAL(prec), INTENT(in)             :: south, north, east, west
   ! Local
   LOGICAL inputProblem

      inputProblem = .FALSE.

      IF( north > 90.0_prec .OR. north < -90.0_prec .OR. &
          south > 90.0_prec .OR. south < -90.0_prec )THEN
         PRINT *, 'Module POP_Regional_Class.f90 : S/R Build :'
         PRINT *, '    Latitudes must be specified in degrees N and'
         PRINT*,  '    must be between -90 and 90.'
         inputProblem = .TRUE.
      ELSE
         myRegion % south = south
         myRegion % north = north
      ENDIF

      IF( east < -360.0_prec .OR. east > 360.0_prec .OR. & 
          west < -360.0_prec .OR. west > 360.0_prec )THEN
         PRINT *, 'Module POP_Regional_Class.f90 : S/R Build :'
         PRINT *, '    Longitudes must be specified in degrees E and'
         PRINT*,  '    must be between -360 and 360.'
         inputProblem = .TRUE.
      ELSE
          
         IF( east < 0.0_prec )THEN
            ! Force the longitudes between 0 and 360 to conform with POP mesh.
            myRegion % east = east + 360.0_prec
         ELSE
            myRegion % east = east
         ENDIF

         IF( west < 0.0_prec )THEN
            ! Force the longitudes between 0 and 360 to conform with POP mesh.
            myRegion % west = west + 360.0_prec
         ELSE
            myRegion % west = west
         ENDIF

      ENDIF

      myRegion % nDOF = mesh % nDOF

      IF( inputProblem )THEN
         STOP '     Stopping!'
      ENDIF
      PRINT*, ' Finding cells in Region '
      CALL myRegion % FindThoseInRegion( mesh )

      PRINT*, ' Finding boundary cells'
      CALL myRegion % FindBoundaryCells( mesh, mystencil, meshType )


 END SUBROUTINE Build_POP_Regional

 SUBROUTINE Trash_POP_Regional( myRegion )
   IMPLICIT NONE
   CLASS( POP_Regional ), INTENT(inout) :: myRegion

      DEALLOCATE( myRegion % ijkInRegion, &
                  myRegion % dofInRegion, &
                  myRegion % inverseDOFMap, &
                  myRegion % boundaryCells )

 END SUBROUTINE Trash_POP_Regional
!
 SUBROUTINE FindThoseInRegion_POP_Regional( myRegion, mesh )

   IMPLICIT NONE
   CLASS( POP_Regional ), INTENT(inout) :: myRegion
   TYPE( POP_Mesh ), INTENT(inout)      :: mesh
   !
   INTEGER :: i, j, k
   INTEGER :: nInRegion
   LOGICAL :: inRegion

      mesh % tracerMask = 1.0_prec
      IF( myRegion % east < myRegion % west  )THEN
         myRegion % crossesPrimeMeridian = .TRUE. 
         PRINT*, '  Region crosses prime-meridian'
         ! In this scenario, the domain crosses over 0 East. As an example
         ! consider the case in the Atlantic Ocean (perhaps near the Agulhas
         ! region), where the western edge of the domain we would like to model
         ! covers longitudes between -10E and 40E. When calling "Build", -10E is
         ! converted to 350, which is greater than 40. To find cells in this
         ! region, we look for cells with lon > 350 and lon < 40. In general,
         ! for this case "lon > mesh % west" OR "lon < mesh % east" 

         ! First, count how many cells are in the region
         nInRegion = 0
         DO j = 1, mesh % nY
            DO i = 1, mesh % nX
               DO k = 1, mesh % KMT(i,j)

                  IF( mesh % tLon(i,j) >= myRegion % west .OR. &
                      mesh % tLon(i,j) <= myRegion % east ) THEN
                     IF( mesh % tLat(i,j) >= myRegion % south .AND. &
                         mesh % tLat(i,j) <= myRegion % north )THEN
                         nInRegion = nInRegion + 1
                     ENDIF
                  ENDIF

               ENDDO
            ENDDO
         ENDDO
         PRINT*, '  Found ', nInRegion, 'in region' 
         myRegion % nCells = nInRegion
         ALLOCATE( myRegion % ijkInRegion(1:3,1:nInRegion), &
                   myRegion % dofInRegion(1:nInRegion), &
                   myRegion % dofToLocalIJK(1:3,nInRegion) )

         nInRegion = 0 
         DO j = 1, mesh % nY
            DO i = 1, mesh % nX
               DO k = 1, mesh % KMT(i,j) 

                  inRegion = .FALSE.
                  IF( mesh % tLon(i,j) >= myRegion % west .OR. &
                      mesh % tLon(i,j) <= myRegion % east ) THEN
                     IF( mesh % tLat(i,j) >= myRegion % south .AND. &
                         mesh % tLat(i,j) <= myRegion % north )THEN

                         nInRegion = nInRegion + 1
                         myRegion % ijkInRegion(1:3,nInRegion) = (/i, j, k/)
                         myRegion % dofInRegion(nInRegion) = mesh % ijkToDOF(i,j,k)

                        inRegion = .TRUE.
                     ENDIF
                  ENDIF

                  IF( .NOT.(inRegion) ) THEN
                     mesh % tracermask(i,j,k) = 0.0_prec
                  ENDIF

               ENDDO
            ENDDO
         ENDDO

      ELSE
         myRegion % crossesPrimeMeridian = .FALSE. 

         ! First, count how many cells are in the region
         nInRegion = 0
         DO j = 1, mesh % nY
            DO i = 1, mesh % nX
               DO k = 1, mesh % KMT(i,j)

                  IF( mesh % tLon(i,j) >= myRegion % west .AND. &
                      mesh % tLon(i,j) <= myRegion % east ) THEN
                     IF( mesh % tLat(i,j) >= myRegion % south .AND. &
                         mesh % tLat(i,j) <= myRegion % north )THEN
                         nInRegion = nInRegion + 1
                     ENDIF
                  ENDIF

               ENDDO
            ENDDO
         ENDDO
         PRINT*, '  Found ', nInRegion, 'in region' 
          
         myRegion % nCells = nInRegion
         ALLOCATE( myRegion % ijkInRegion(1:3,1:nInRegion), &
                   myRegion % dofInRegion(1:nInRegion), &
                   myRegion % doftoLocalIJK(1:3,1:nInRegion) )

         nInRegion = 0 
         DO j = 1, mesh % nY
            DO i = 1, mesh % nX
               DO k = 1, mesh % KMT(i,j) 

                  inRegion = .FALSE.
                  IF( mesh % tLon(i,j) >= myRegion % west .AND. &
                      mesh % tLon(i,j) <= myRegion % east ) THEN
                     IF( mesh % tLat(i,j) >= myRegion % south .AND. &
                         mesh % tLat(i,j) <= myRegion % north )THEN

                         nInRegion = nInRegion + 1
                         myRegion % ijkInRegion(1:3,nInRegion) = (/i, j, k/)
                         myRegion % dofInRegion(nInRegion) = mesh % ijkToDOF(i,j,k)
                         inRegion = .TRUE.

                     ENDIF
                  ENDIF

                  IF( .NOT.(inRegion) ) THEN
                     mesh % tracermask(i,j,k) = 0.0_prec
                  ENDIF

               ENDDO
            ENDDO
         ENDDO

      ENDIF

      ! Generate the inverse dof map
      ALLOCATE( myRegion % inverseDOFMap(1:mesh % nDOF) )
      myRegion % inverseDOFMap = 0

      ! Determine the inverse map
      DO k = 1, myRegion % nCells
         ! k is the local DOF, dofInRegion(k) is the globalDOF
         myRegion % inverseDOFMap( myRegion % dofInRegion(k) ) = k
      ENDDO

                  
 END SUBROUTINE FindThoseInRegion_POP_Regional
!
 SUBROUTINE FindBoundaryCells_POP_Regional( myRegion, mesh, mystencil, meshType )

   IMPLICIT NONE
   CLASS( POP_Regional ), INTENT(inout)  :: myRegion
   TYPE( POP_Mesh ), INTENT(inout)       :: mesh
   TYPE( Stencil ), INTENT(in)           :: mystencil
   INTEGER, INTENT(in)                   :: meshType
   ! Local
   INTEGER :: i, j, k, ii, jj, this_i, this_j, m, n
   INTEGER :: minI, maxI, minJ, maxJ, minL
   INTEGER :: nBCells
   INTEGER :: regionmask(1:mesh % nX,1:mesh % nY,1:mesh % nZ)
   INTEGER :: stencilSum
      
      ! To find the boundary cells, we will create a "mask" field with 0's assigned to the 
      ! points within the specified region and one assigned to points outside.
      ! Then, for every point in the region, we add the values of the neighbors that lie within the
      ! stencil.
      
      ! (1) Generate the regional mask
      regionMask = 1 ! Set all values to 1 initially

      minI = MINVAL( myRegion % ijkInRegion(1,:) ) 
      maxI = MAXVAL( myRegion % ijkInRegion(1,:) ) 
      minJ = MINVAL( myRegion % ijkInRegion(2,:) ) 
      maxJ = MAXVAL( myRegion % ijkInRegion(2,:) ) 

      IF( myRegion % crossesPrimeMeridian .AND. mesh % tLon(minI,minJ) > myRegion % west )THEN
         minL = mesh % tLon(minI,minJ) - 360.0_prec
      ELSE
         minL = mesh % tLon(minI,minJ)
      ENDIF

      ! Zeros out the mask points within the region
      IF( mesh % tLon(maxI,minJ) < minL )THEN
         
         DO k = 1, mesh % nZ
            DO j = minJ, maxJ
               DO i = 1, minI
                  regionMask(i,j,k) = 0
               ENDDO
            ENDDO
         ENDDO

         DO k = 1, mesh % nZ
            DO j = minJ, maxJ
               DO i = maxI, mesh % nX
                  regionMask(i,j,k) = 0
               ENDDO
            ENDDO
         ENDDO

      ELSE

         DO k = 1, mesh % nZ
            DO j = minJ, maxJ
               DO i = minI, maxI
                  regionMask(i,j,k) = 0
               ENDDO
            ENDDO
         ENDDO

      ENDIF

      ! (2) Loop over the region and compute the sum of the points within the stencil. The first 
      ! pass through will count the number of border cells.
      nBCells = 0
      DO m = 1, myRegion % nCells

         i = myRegion % ijkInRegion(1,m) ! Global i, j, k
         j = myRegion % ijkInRegion(2,m)
         k = myRegion % ijkInRegion(3,m)

         stencilSum = 0
         DO n = 1, myStencil % nPoints

            ii = i + myStencil % relativeNeighbors(1,n)
            jj = j + myStencil % relativeNeighbors(2,n)
            CALL GetTrueIJ( meshType, ii, jj, &
                            mesh % nX, mesh % nY, &
                            this_i, this_j )
           
            stencilSum = stencilSum + regionMask(this_i,this_j,k)

         ENDDO

         IF( stencilSum > 0 )THEN
            nBCells = nBcells + 1
            mesh % tracerMask(i,j,k) = -1.0_prec ! Set the border cell tracer mask to -1
         ENDIF

      ENDDO

      PRINT*,'  Found', nBCells,' boundary cells.' 
      myRegion % nBCells = nBCells
      ALLOCATE( myRegion % boundaryCells(1:nBCells) )
      ! (3) Use the mesh tracermask to fill in the borderCells attribute of the POP_Regional data structure
      nBCells = 0
      DO m = 1, myRegion % nCells

         i = myRegion % ijkInRegion(1,m)
         j = myRegion % ijkInRegion(2,m)
         k = myRegion % ijkInRegion(3,m)

         IF( mesh % tracerMask(i,j,k) == -1.0_prec )THEN
            nBCells = nBCells + 1
            myRegion % boundaryCells(nBCells) = m
         ENDIF

      ENDDO

 END SUBROUTINE FindBoundaryCells_POP_Regional
!
 SUBROUTINE GenerateRegionalMesh_POP_Regional( myRegion, mesh, regionalmesh )

   IMPLICIT NONE
   CLASS( POP_Regional ), INTENT(inout) :: myRegion
   TYPE( POP_Mesh ), INTENT(in)         :: mesh
   TYPE( POP_Mesh ), INTENT(out)        :: regionalMesh
   ! Local
   INTEGER :: minI, maxI, minJ, maxJ, i, j, k, m, ii, jj, nXr, nYr
   REAL(prec) :: minL

      minI = MINVAL( myRegion % ijkInRegion(1,:) ) 
      maxI = MAXVAL( myRegion % ijkInRegion(1,:) ) 
      minJ = MINVAL( myRegion % ijkInRegion(2,:) ) 
      maxJ = MAXVAL( myRegion % ijkInRegion(2,:) ) 

      IF( myRegion % crossesPrimeMeridian .AND. mesh % tLon(minI,minJ) > myRegion % west )THEN
         minL = mesh % tLon(minI,minJ) - 360.0_prec
      ELSE
         minL = mesh % tLon(minI,minJ)
      ENDIF


      IF( mesh % tLon(maxI,minJ) < minL )THEN
         ! In this case the region crosses the boundary of the periodic computational mesh
         nXr = ( mesh % nX - maxI + 1 ) + minI
         nYr = maxJ - minJ + 1

         CALL regionalMesh % Build( nXr, nYr, mesh % nZ )

        m = 0
        jj = 0
        DO j = minJ, maxJ
           jj = jj+1
           ii = 0
           DO i = maxI, mesh % nX
              ii = ii + 1
              regionalMesh % tLon(ii,jj)  = mesh % tLon(i,j)
              regionalMesh % tLat(ii,jj)  = mesh % tLat(i,j)
              regionalMesh % dXt(ii,jj)   = mesh % dXt(i,j)
              regionalMesh % dYt(ii,jj)   = mesh % dYt(i,j)
              regionalMesh % tArea(ii,jj) = mesh % tArea(i,j)
              regionalMesh % KMT(ii,jj)   = mesh % KMT(i,j)
              DO k = 1, mesh % KMT(i,j)
              regionalMesh % tracerMask(ii,jj,k)  = mesh % tracerMask(i,j,k)
              IF( mesh % tracerMask(i,j,k) /= 0.0_prec )THEN
                 m = m+1
                 myRegion % dofToLocalIJK(1:3,m) = (/ ii, jj, k /)
              ENDIF
              ENDDO
           ENDDO

           DO i = 1, minI
              ii = ii + 1
              regionalMesh % tLon(ii,jj)  = mesh % tLon(i,j)
              regionalMesh % tLat(ii,jj)  = mesh % tLat(i,j)
              regionalMesh % dXt(ii,jj)   = mesh % dXt(i,j)
              regionalMesh % dYt(ii,jj)   = mesh % dYt(i,j)
              regionalMesh % tArea(ii,jj) = mesh % tArea(i,j)
              regionalMesh % KMT(ii,jj)   = mesh % KMT(i,j)
              DO k = 1, mesh % KMT(i,j)
              regionalMesh % tracerMask(ii,jj,k)  = mesh % tracerMask(i,j,k)
              IF( mesh % tracerMask(i,j,k) /= 0.0_prec )THEN
                 m = m+1
                 myRegion % dofToLocalIJK(1:3,m) = (/ ii,jj,k /)
              ENDIF
              ENDDO
           ENDDO

        ENDDO

     ELSE

         nXr = maxI - minI + 1
         nYr = maxJ - minJ + 1

         CALL regionalMesh % Build( nXr, nYr, mesh % nZ )

        jj = 0
        DO j = minJ, maxJ
           jj = jj+1
           ii = 0
           DO i = minI, maxI
              ii = ii + 1
              regionalMesh % tLon(ii,jj)  = mesh % tLon(i,j)
              regionalMesh % tLat(ii,jj)  = mesh % tLat(i,j)
              regionalMesh % dXt(ii,jj)   = mesh % dXt(i,j)
              regionalMesh % dYt(ii,jj)   = mesh % dYt(i,j)
              regionalMesh % tArea(ii,jj) = mesh % tArea(i,j)
              regionalMesh % KMT(ii,jj)   = mesh % KMT(i,j)
              DO k = 1, mesh % KMT(i,j)
              regionalMesh % tracerMask(ii,jj,k)  = mesh % tracerMask(i,j,k)
              IF( mesh % tracerMask(i,j,k) /= 0.0_prec )THEN
                 m = m+1
                 myRegion % dofToLocalIJK(1:3,m) = (/ ii,jj,k /)
              ENDIF
              ENDDO
           ENDDO

        ENDDO

     ENDIF

     IF( myRegion % crossesPrimeMeridian )THEN

        DO j = 1, nYr
           DO i = 1, nXr
              IF( regionalMesh % tLon(i,j) >= myRegion % west )THEN
                 regionalMesh % tLon(i,j) = regionalMesh % tLon(i,j) - 360.0_prec
              ENDIF
           ENDDO
        ENDDO

     ENDIF

     DO k = 1, mesh % nZ
        regionalMesh % z(k)   = mesh % z(k)
        regionalMesh % dz(k)  = mesh % dz(k)
        regionalMesh % dzw(k) = mesh % dzw(k)
     ENDDO


 END SUBROUTINE GenerateRegionalMesh_POP_Regional
!
 SUBROUTINE WritePickup_POP_Regional( myRegional, filename )

   IMPLICIT NONE
   CLASS( POP_Regional ), INTENT(in) :: myRegional
   CHARACTER(*), INTENT(in)          :: filename
   ! Local 
   INTEGER :: fUnit, i


      OPEN( UNIT=NewUnit(fUnit), &
            FILE=TRIM(filename)//'.regional',&
            FORM='FORMATTED',&
            ACCESS='SEQUENTIAL',&
            STATUS='REPLACE', &
            ACTION='WRITE' )

     WRITE( fUnit, * ) myRegional % nCells, myRegional % nBCells, myRegional % nDOF
     
     DO i = 1, myRegional % nCells
        WRITE(fUnit,*) myRegional % ijkInRegion(1:3,i), &
                       myRegional % dofInRegion(i), &
                       myRegional % dofToLocalIJK(1:3,i)
     ENDDO

     DO i = 1, myRegional % nBCells
        WRITE(fUnit,*)myRegional % boundaryCells(i)
     ENDDO

     DO i = 1, myRegional % nDOF
        WRITE(fUnit,*)myRegional % inverseDOFMap(i)
     ENDDO

 END SUBROUTINE WritePickup_POP_Regional
!
 SUBROUTINE ReadPickup_POP_Regional( myRegional, filename )

   IMPLICIT NONE
   CLASS( POP_Regional ), INTENT(inout) :: myRegional
   CHARACTER(*), INTENT(in)             :: filename
   ! Local 
   INTEGER :: fUnit, i, readStatus


      OPEN( UNIT=NewUnit(fUnit), &
            FILE=TRIM(filename)//'.regional',&
            FORM='FORMATTED',&
            ACCESS='SEQUENTIAL',&
            STATUS='OLD', &
            ACTION='READ' )

     READ( fUnit, * ) myRegional % nCells, myRegional % nBCells, myRegional % nDOF
     
     ALLOCATE( myRegional % ijkInRegion(1:3,1:myRegional % nCells), &
               myRegional % dofInRegion(1:myRegional % nCells), &
               myRegional % dofToLocalIJK(1:3,1:myRegional % nCells), &
               myRegional % boundaryCells(1:myRegional % nBCells), &
               myRegional % inverseDOFMap(1:myRegional % nDOF) )

     DO i = 1, myRegional % nCells
        READ(fUnit,*) myRegional % ijkInRegion(1:3,i), &
                       myRegional % dofInRegion(i), &
                       myRegional % dofToLocalIJK(1:3,i)
     ENDDO

     DO i = 1, myRegional % nBCells
        READ(fUnit,*)myRegional % boundaryCells(i)
     ENDDO

     DO i = 1, myRegional % nDOF
        READ(fUnit,*,IOSTAT=readStatus)myRegional % inverseDOFMap(i)
        IF(readStatus < 0 )THEN
           PRINT*, 'End of File at DOF : ',i
        ENDIF
     ENDDO

 END SUBROUTINE ReadPickup_POP_Regional


END MODULE POP_Regional_Class
