!  POP_Regional_Class.f90
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


MODULE POP_Regional_Class


! src/common/
USE ModelPrecision
! src/POP/
USE POP_Mesh_Class
USE POP_Stencil_Class
USE POP_GridTypeMappings
!
USE netcdf

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

   TYPE BoundaryMap
      INTEGER :: nBCells, nPCells
      INTEGER, ALLOCATABLE :: boundaryCells(:)
      INTEGER, ALLOCATABLE :: prescribedCells(:)
      
   END TYPE BoundaryMap

   TYPE POP_Regional
      REAL(prec) :: south, north, east, west
      LOGICAL    :: crossesPrimeMeridian
      INTEGER    :: nCells, nDOF, nMasks
      INTEGER, ALLOCATABLE :: ijkInRegion(:,:)
      INTEGER, ALLOCATABLE :: dofInRegion(:)  
      INTEGER, ALLOCATABLE :: inverseDOFMap(:)  
      INTEGER, ALLOCATABLE :: dofToLocalIJK(:,:)
      TYPE( BoundaryMap ), ALLOCATABLE  :: bMap(:)

      CONTAINS

      PROCEDURE :: Build => Build_POP_Regional
      PROCEDURE :: Trash => Trash_POP_Regional
      PROCEDURE :: LoadMaskField => LoadMaskField_POP_Regional

      PROCEDURE :: FindThoseInRegion => FindThoseInRegion_POP_Regional 
      PROCEDURE :: FindBoundaryCells => FindBoundaryCells_POP_Regional

      PROCEDURE :: GenerateRegionalMesh => GenerateRegionalMesh_POP_Regional

      PROCEDURE :: WritePickup => WritePickup_POP_Regional
      PROCEDURE :: ReadPickup  => ReadPickup_POP_Regional

   END TYPE POP_Regional


CONTAINS


 SUBROUTINE Build_POP_Regional( myRegion, mesh, regionalMesh, mystencil, meshType, south, north, east, west, maskfile )

   IMPLICIT NONE
   CLASS( POP_Regional ), INTENT(out) :: myRegion
   TYPE( POP_Mesh ), INTENT(inout)    :: mesh
   TYPE( POP_Mesh ), INTENT(inout)    :: regionalMesh
   TYPE( Stencil ), INTENT(in)        :: mystencil
   INTEGER, INTENT(in)                :: meshType
   REAL(prec), INTENT(in)             :: south, north, east, west
   CHARACTER(*), INTENT(in)           :: maskfile
   !INTEGER, INTENT(in), OPTIONAL      :: maskfield(1:mesh % nX,1:mesh % nY,1:myRegion % nMasks)
   ! Local
   LOGICAL              :: inputProblem
   INTEGER, ALLOCATABLE :: maskfield(:,:,:)

      inputProblem = .FALSE.
      IF( TRIM(maskfile) /= '' )THEN
      !IF( present(maskfield) )THEN

         CALL myRegion % LoadMaskField( mesh, maskfield, maskfile ) 
         ALLOCATE( myRegion % bMap(1:myRegion % nMasks) )

         PRINT*, ' Finding cells in Region '
         CALL myRegion % FindThoseInRegion( mesh, maskfield )

         PRINT*, ' Finding boundary cells'
         CALL myRegion % FindBoundaryCells( mesh, mystencil, meshType, maskfield )

         CALL myRegion % GenerateRegionalMesh( mesh, regionalMesh, maskfield )
         DEALLOCATE( maskfield )

      ELSE
         myRegion % nMasks = 1
         ALLOCATE( myRegion % bMap(1:myRegion % nMasks) )

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

         CALL myRegion % GenerateRegionalMesh( mesh, regionalMesh )

      ENDIF

 END SUBROUTINE Build_POP_Regional

 SUBROUTINE Trash_POP_Regional( myRegion )
   IMPLICIT NONE
   CLASS( POP_Regional ), INTENT(inout) :: myRegion
   ! Local
   INTEGER :: i

      DEALLOCATE( myRegion % ijkInRegion, &
                  myRegion % dofInRegion, &
                  myRegion % inverseDOFMap )

      DO i = 1, myRegion % nMasks
         DEALLOCATE( myRegion % bMap(i) % boundaryCells, &
                     myRegion % bMap(i) % prescribedCells )
      ENDDO
      DEALLOCATE( myRegion % bMap )

 END SUBROUTINE Trash_POP_Regional
!
 SUBROUTINE LoadMaskField_POP_Regional( myRegion, mesh, maskfield, maskfile )

   IMPLICIT NONE
   CLASS( POP_Regional ), INTENT(inout) :: myRegion
   TYPE( POP_Mesh ), INTENT(inout)      :: mesh
   INTEGER, ALLOCATABLE, INTENT(out)    :: maskfield(:,:,:)
   CHARACTER(*), INTENT(in)             :: maskfile
   ! Local
   INTEGER      :: start(1:2), recCount(1:2)
   INTEGER      :: ncid, varid, i
   CHARACTER(3) :: maskChar

      start    = (/1, 1/)
      recCount = (/mesh % nX, mesh % nY/)

      CALL Check( nf90_open( TRIM(maskfile), nf90_nowrite, ncid ) )

      CALL Check( nf90_inq_varid( ncid, "nMasks", varid ) )
      CALL Check( nf90_get_var( ncid,  varid, myRegion % nMasks ) )

      ALLOCATE( maskfield(1:mesh % nX, 1:mesh % nY, 1:myRegion % nMasks) )

      DO i = 1, myRegion % nMasks

         WRITE( maskChar, '(I3.3)')i
         CALL Check( nf90_inq_varid( ncid, "mask"//maskChar, varid ) )
         CALL Check( nf90_get_var( ncid,  varid, maskfield(:,:,i), start, recCount ) )

      ENDDO

      CALL Check( nf90_close( ncid ) )


 END SUBROUTINE LoadMaskField_POP_Regional
!
 SUBROUTINE FindThoseInRegion_POP_Regional( myRegion, mesh, maskfield )

   IMPLICIT NONE
   CLASS( POP_Regional ), INTENT(inout) :: myRegion
   TYPE( POP_Mesh ), INTENT(inout)      :: mesh
   INTEGER, INTENT(in), OPTIONAL        :: maskfield(1:mesh % nX, 1:mesh % nY,1:myRegion % nMasks)
   !
   INTEGER :: i, j, k
   INTEGER :: nInRegion
   LOGICAL :: inRegion
   LOGICAL :: logicMask(1:mesh % nX, 1:mesh % nY)
   REAL(prec) :: minLonMesh, maxLonMesh, minLonRegion, maxLonRegion

      IF( present(maskfield) )THEN
        
          logicMask = .FALSE.

          mesh % tracerMask = fillValue

            nInRegion = 0
            DO j = 1, mesh % nY
               DO i = 1, mesh % nX
                  DO k = 1, mesh % KMT(i,j)
   
                     IF( maskfield(i,j,1) /= 0 )THEN
                         nInRegion = nInRegion + 1
                         mesh % tracerMask(i,j,k) = 1.0_prec
                         logicMask(i,j) = .TRUE.
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
   
                     IF( maskfield(i,j,1) /= 0 )THEN
                         nInRegion = nInRegion + 1
                         myRegion % ijkInRegion(1:3,nInRegion) = (/i, j, k/)
                         myRegion % dofInRegion(nInRegion) = mesh % ijkToDOF(i,j,k)
                     ENDIF
   
                  ENDDO
               ENDDO
            ENDDO
         ! Generate the inverse dof map
         ALLOCATE( myRegion % inverseDOFMap(1:mesh % nDOF) )
         myRegion % inverseDOFMap = 0
   
         ! Determine the inverse map
         DO k = 1, myRegion % nCells
            ! k is the local DOF, dofInRegion(k) is the globalDOF
            myRegion % inverseDOFMap( myRegion % dofInRegion(k) ) = k
         ENDDO
         
         myRegion % south = MINVAL( mesh % tLat, logicMask )
         myRegion % north = MAXVAL( mesh % tLat, logicMask )

         ! Need to determine if the region crosses the prime meridian.
         ! To do this, we need to find the min and max longitudes. If the
         ! minimum longitude equals the minimum longitude in the mesh and the
         ! maximum longitude and the maximum longitude equals the maximum
         ! longitude in the mesh, then the domain crosses the prime meridian.
         minLonMesh = MINVAL( mesh % tLon ) 
         maxLonMesh = MAXVAL( mesh % tLon ) 
         minLonRegion = MINVAL( mesh % tLon, logicMask )
         maxLonRegion = MAXVAL( mesh % tLon, logicMask )

         IF( minLonRegion == minLonMesh .AND. maxLonRegion == maxLonMesh )THEN

            PRINT*, '  Region crosses prime-meridian'
            myRegion % crossesPrimeMeridian = .TRUE.

            ! Shift the mesh
            DO j = 1, mesh % nY
               DO i = 1, mesh % nX
                  IF( mesh % tLon(i,j) >= 180.0_prec )THEN
                     mesh % tLon(i,j) = mesh % tLon(i,j) - 360.0_prec
                  ENDIF
               ENDDO
           ENDDO

           minLonMesh = MINVAL( mesh % tLon ) 
           maxLonMesh = MAXVAL( mesh % tLon ) 
           minLonRegion = MINVAL( mesh % tLon, logicMask )
           maxLonRegion = MAXVAL( mesh % tLon, logicMask )
 
           IF( minLonRegion == minLonMesh .AND. maxLonRegion == maxLonMesh )THEN
              ! In this case, the region completely encompasses the globe in a
              ! latitudinal band.
              PRINT*, '  Wrap-around detected'
              DO j = 1, mesh % nY
                 DO i = 1, mesh % nX
                    IF( mesh % tLon(i,j) < 0.0_prec )THEN
                       mesh % tLon(i,j) = mesh % tLon(i,j) + 360.0_prec
                    ENDIF
                 ENDDO
             ENDDO
             myRegion % west = minLonMesh
             myRegion % east = maxLonMesh
           ELSE

             myRegion % west = minLonRegion + 360.0_prec
             myRegion % east = maxLonRegion

           ENDIF

        ELSE

           myRegion % crossesPrimeMeridian = .TRUE.
           myRegion % west = minLonRegion
           myRegion % east = maxLonRegion
        ENDIF

      ELSE
   
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

      ENDIF
         
 END SUBROUTINE FindThoseInRegion_POP_Regional
!
 SUBROUTINE FindBoundaryCells_POP_Regional( myRegion, mesh, mystencil, meshType, maskfield )

   IMPLICIT NONE
   CLASS( POP_Regional ), INTENT(inout)  :: myRegion
   TYPE( POP_Mesh ), INTENT(inout)       :: mesh
   TYPE( Stencil ), INTENT(in)           :: mystencil
   INTEGER, INTENT(in)                   :: meshType
   INTEGER, INTENT(in), OPTIONAL         :: maskfield(1:mesh % nX, 1:mesh % nY,1:myRegion % nMasks)
   ! Local
   INTEGER :: i, j, k, ii, jj, this_i, this_j, m, n, iMask
   INTEGER :: minI, maxI, minJ, maxJ
   REAL(prec) :: minL
   INTEGER :: nBCells, nPCells
   INTEGER :: regionmask(1:mesh % nX,1:mesh % nY,1:mesh % nZ,1:myRegion % nMasks)
   INTEGER :: stencilSum
      
      IF( present(maskfield) )THEN

         ! To find the boundary cells, we will create a "mask" field with 0's assigned to the 
         ! points within the specified region and one assigned to points outside.
         ! Then, for every point in the region, we add the values of the neighbors that lie within the
         ! stencil.
         
         ! (1) Generate the regional mask
         regionMask = 1 ! Set all values to 1 initially
         DO iMask = 1, myRegion % nMasks

         DO k = 1, mesh % nZ
            DO j = 1, mesh % nY
               DO i = 1, mesh % nX
                  IF( maskfield(i,j,iMask) /= 0 )THEN
                     regionMask(i,j,k,iMask) = 0
                  ENDIF
               ENDDO
            ENDDO
         ENDDO

         ! (2) Loop over the region and compute the sum of the points within the stencil. The first 
         ! pass through will count the number of border cells.
         nBCells = 0
         nPCells = 0
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
              
               stencilSum = stencilSum + regionMask(this_i,this_j,k,iMask)

            ENDDO

            IF( stencilSum > 0 .OR. maskfield(i,j,iMask) == -1 )THEN
               nBCells = nBcells + 1
               IF( maskfield(i,j,iMask) == -1 )THEN ! In this case the user wants to prescribe a boundary condition
                  nPCells = nPCells + 1
                  mesh % tracerMask(i,j,k) = -REAL(iMask,prec)
               ELSE
                  IF( mesh % tracerMask(i,j,k) > 0.0_prec )THEN
                     mesh % tracerMask(i,j,k) = 0.0_prec
                  ENDIF
               ENDIF
            ENDIF

         ENDDO

         PRINT*,'  Found', nBCells,' boundary cells for Mask ID',iMask,'.' 
         PRINT*,'  Found', nPCells,' prescribed cells for Mask ID',iMask,'.' 
         myRegion % bMap(iMask) % nBCells = nBCells
         myRegion % bMap(iMask) % nPCells = nPCells
         ALLOCATE( myRegion % bMap(iMask) % boundaryCells(1:nBCells) )
         ALLOCATE( myRegion % bMap(iMask) % prescribedCells(1:nBCells) )
         ! (3) Use the mesh tracermask to fill in the borderCells attribute of the POP_Regional data structure
         nBCells = 0
         nPCells = 0
         DO m = 1, myRegion % nCells

            i = myRegion % ijkInRegion(1,m)
            j = myRegion % ijkInRegion(2,m)
            k = myRegion % ijkInRegion(3,m)

            IF( mesh % tracerMask(i,j,k) == 0.0_prec )THEN
               nBCells = nBCells + 1
               myRegion % bMap(iMask) % boundaryCells(nBCells) = m
            ELSEIF( mesh % tracerMask(i,j,k) == -REAL(iMask,prec) )THEN
               nPCells = nPCells + 1
               nBCells = nBCells + 1
               myRegion % bMap(iMask) % boundaryCells(nBCells)   = m
               myRegion % bMap(iMask) % prescribedCells(nPCells) = m
            ENDIF

         ENDDO

         ENDDO

      ELSE
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
                     regionMask(i,j,k,1) = 0
                  ENDDO
               ENDDO
            ENDDO

            DO k = 1, mesh % nZ
               DO j = minJ, maxJ
                  DO i = maxI, mesh % nX
                     regionMask(i,j,k,1) = 0
                  ENDDO
               ENDDO
            ENDDO

         ELSE

            DO k = 1, mesh % nZ
               DO j = minJ, maxJ
                  DO i = minI, maxI
                     regionMask(i,j,k,1) = 0
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
              
               stencilSum = stencilSum + regionMask(this_i,this_j,k,1)

            ENDDO

            IF( stencilSum > 0 )THEN
               nBCells = nBcells + 1
               mesh % tracerMask(i,j,k) = -1.0_prec ! Set the border cell tracer mask to -1
            ENDIF

         ENDDO

         PRINT*,'  Found', nBCells,' boundary cells.' 
         myRegion % bMap(1) %  nBCells = nBCells
         ALLOCATE( myRegion % bMap(1) % boundaryCells(1:nBCells) )
         ALLOCATE( myRegion % bMap(1) % prescribedCells(1) )
         ! (3) Use the mesh tracermask to fill in the borderCells attribute of the POP_Regional data structure
         nBCells = 0
         DO m = 1, myRegion % nCells

            i = myRegion % ijkInRegion(1,m)
            j = myRegion % ijkInRegion(2,m)
            k = myRegion % ijkInRegion(3,m)

            IF( mesh % tracerMask(i,j,k) == -1.0_prec )THEN
               nBCells = nBCells + 1
               myRegion % bMap(1) % boundaryCells(nBCells) = m
            ENDIF

         ENDDO

      ENDIF

 END SUBROUTINE FindBoundaryCells_POP_Regional
!
 SUBROUTINE GenerateRegionalMesh_POP_Regional( myRegion, mesh, regionalmesh, maskfield )

   IMPLICIT NONE
   CLASS( POP_Regional ), INTENT(inout) :: myRegion
   TYPE( POP_Mesh ), INTENT(in)         :: mesh
   TYPE( POP_Mesh ), INTENT(out)        :: regionalMesh
   INTEGER, INTENT(in), OPTIONAL        :: maskfield(1:mesh % nX, 1:mesh % nY,1:myRegion % nMasks)
   ! Local
   INTEGER :: minI, maxI, minJ, maxJ, i, j, k, m, ii, jj, nXr, nYr
   REAL(prec) :: minL

      IF( present(maskfield) )THEN
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
                 IF( maskfield(i,j,1) /= 0 )THEN
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
                 IF( maskfield(i,j,1) /= 0 )THEN
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
           m  = 0   
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
   
      ELSE
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
   
           m = 0
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
     ENDIF


 END SUBROUTINE GenerateRegionalMesh_POP_Regional
!
 SUBROUTINE WritePickup_POP_Regional( myRegional, filename, maskProvided )

   IMPLICIT NONE
   CLASS( POP_Regional ), INTENT(in) :: myRegional
   CHARACTER(*), INTENT(in)          :: filename
   LOGICAL, INTENT(in)               :: maskProvided
   ! Local 
   INTEGER :: fUnit, i, j


      OPEN( UNIT=NewUnit(fUnit), &
            FILE=TRIM(filename)//'.regional',&
            FORM='FORMATTED',&
            ACCESS='SEQUENTIAL',&
            STATUS='REPLACE', &
            ACTION='WRITE' )

     IF( maskProvided )THEN
        WRITE( fUnit, * ) myRegional % nCells, myRegional % nMasks, myRegional % nDOF
     ELSE
        WRITE( fUnit, * ) myRegional % nCells, myRegional % bMap(1) % nBCells, myRegional % nDOF
     ENDIF    
 
     DO i = 1, myRegional % nCells
        WRITE(fUnit,*) myRegional % ijkInRegion(1:3,i), &
                       myRegional % dofInRegion(i), &
                       myRegional % dofToLocalIJK(1:3,i)
     ENDDO

     DO j = 1, myRegional % nMasks
        WRITE( fUnit, * ) myRegional % bMap(j) % nBCells
        DO i = 1, myRegional % bMap(j) % nBCells
           WRITE(fUnit,*)myRegional % bMap(j) % boundaryCells(i)
        ENDDO
     ENDDO

     IF( maskProvided )THEN
        DO j = 1, myRegional % nMasks
           WRITE( fUnit, * ) myRegional % bMap(j) % nPCells
           DO i = 1, myRegional % bMap(j) % nPCells
              WRITE(fUnit,*)myRegional % bMap(j) % prescribedCells(i)
           ENDDO
        ENDDO
     ENDIF

     DO i = 1, myRegional % nDOF
        WRITE(fUnit,*)myRegional % inverseDOFMap(i)
     ENDDO

 END SUBROUTINE WritePickup_POP_Regional
!
 SUBROUTINE ReadPickup_POP_Regional( myRegional, filename, maskProvided )

   IMPLICIT NONE
   CLASS( POP_Regional ), INTENT(inout) :: myRegional
   CHARACTER(*), INTENT(in)             :: filename
   LOGICAL, INTENT(in)                  :: maskProvided
   ! Local 
   INTEGER :: fUnit, i, j, readStatus


      OPEN( UNIT=NewUnit(fUnit), &
            FILE=TRIM(filename)//'.regional',&
            FORM='FORMATTED',&
            ACCESS='SEQUENTIAL',&
            STATUS='OLD', &
            ACTION='READ' )

     IF( maskProvided )THEN
        READ( fUnit, * ) myRegional % nCells, myRegional % nMasks, myRegional % nDOF
     
        ALLOCATE( myRegional % ijkInRegion(1:3,1:myRegional % nCells), &
                  myRegional % dofInRegion(1:myRegional % nCells), &
                  myRegional % dofToLocalIJK(1:3,1:myRegional % nCells), &
                  myRegional % bMap(1:myRegional % nMasks), &
                  myRegional % inverseDOFMap(1:myRegional % nDOF) )
     ELSE

        ALLOCATE( myRegional % bMap(1) )
        READ( fUnit, * ) myRegional % nCells, myRegional % bMap(1) % nBCells, myRegional % nDOF
        ALLOCATE( myRegional % ijkInRegion(1:3,1:myRegional % nCells), &
                  myRegional % dofInRegion(1:myRegional % nCells), &
                  myRegional % dofToLocalIJK(1:3,1:myRegional % nCells), &
                  myRegional % inverseDOFMap(1:myRegional % nDOF) )
        myRegional %  bMap(1) % nPCells = 0
        myRegional % nMasks  = 1
        ALLOCATE( myRegional % bMap(1) % boundaryCells(1:myRegional % bMap(1) % nBCells),&
                  myRegional % bMap(1) % prescribedCells(1) )
     ENDIF    
     
     DO i = 1, myRegional % nCells
        READ(fUnit,*) myRegional % ijkInRegion(1:3,i), &
                       myRegional % dofInRegion(i), &
                       myRegional % dofToLocalIJK(1:3,i)
     ENDDO

     DO j = 1, myRegional % nMasks
        READ( fUnit, * ) myRegional % bMap(j) % nBCells
 
        IF( .NOT. ALLOCATED( myRegional % bMap(j) % boundaryCells ) )THEN
          ALLOCATE( myRegional % bMap(j) % boundaryCells(1:myRegional % bMap(j) % nBCells) )
        ENDIF

        DO i = 1, myRegional % bMap(j) % nBCells
           READ(fUnit,*)myRegional % bMap(j) % boundaryCells(i)
        ENDDO

     ENDDO

     IF( maskProvided )THEN
        DO j = 1, myRegional % nMasks
           READ( fUnit, * ) myRegional % bMap(j) % nPCells
           ALLOCATE( myRegional % bMap(j) % prescribedCells(1:myRegional % bMap(j) % nPCells) )
           DO i = 1, myRegional % bMap(j) % nPCells
              READ(fUnit,*)myRegional % bMap(j) % prescribedCells(i)
           ENDDO
        ENDDO
     ENDIF

     DO i = 1, myRegional % nDOF
        READ(fUnit,*,IOSTAT=readStatus)myRegional % inverseDOFMap(i)
        IF(readStatus < 0 )THEN
           PRINT*, 'End of File at DOF : ',i
        ENDIF
     ENDDO

 END SUBROUTINE ReadPickup_POP_Regional


END MODULE POP_Regional_Class
