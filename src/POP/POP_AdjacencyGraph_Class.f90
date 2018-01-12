! POP_AdjacencyGraph_Class.f90
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


MODULE POP_AdjacencyGraph_Class


! src/common/
USE ModelPrecision
USE CommonRoutines
! src/POP/
USE POP_Stencil_Class
USE POP_Mesh_Class
USE POP_GridTypeMappings


IMPLICIT NONE

   TYPE POP_AdjacencyGraph
      INTEGER(KIND=8)      :: nDOF
      INTEGER              :: maxValence
      INTEGER              :: nColors
      INTEGER, ALLOCATABLE :: valence(:)
      INTEGER, ALLOCATABLE :: color(:)
      INTEGER(KIND=8), ALLOCATABLE :: neighbors(:,:)

      CONTAINS

      PROCEDURE :: Build => Build_POP_AdjacencyGraph      
      PROCEDURE :: Trash => Trash_POP_AdjacencyGraph      

      PROCEDURE :: ConstructFromStencil => ConstructFromStencil_POP_AdjacencyGraph
      PROCEDURE :: GreedyColoring       => GreedyColoring_POP_AdjacencyGraph
    
      PROCEDURE :: ReadGraphFile     => ReadGraphFile_POP_AdjacencyGraph
      PROCEDURE :: WriteGraphFile    => WriteGraphFile_POP_AdjacencyGraph
      PROCEDURE :: ReadGraphBinFile  => ReadGraphBinFile_POP_AdjacencyGraph
      PROCEDURE :: WriteGraphBinFile => WriteGraphBinFile_POP_AdjacencyGraph


   END TYPE POP_AdjacencyGraph


CONTAINS

 SUBROUTINE Build_POP_AdjacencyGraph( myGraph, nDOF, maxValence )
 ! This subroutine allocates space for the adjacency graph and initializes
 ! all attributes to 0.
 ! INPUT: 
 !     nDOF - is the number of Degrees Of Freedom; ie, the number of wet
 !     gridpoints in the POP mesh.
 !     
 !     maxValence - the maximum number of neighbors that a grid point will
 !     have. This number depends on the stencil that is being used to construct
 !     the adjacency graph.
 ! ============================================================================ !
   IMPLICIT NONE
   CLASS( POP_AdjacencyGraph ), INTENT(inout) :: myGraph
   INTEGER, INTENT(in)                        :: nDOF
   INTEGER, INTENT(in)                        :: maxValence
   
      myGraph % nDOF       = nDOF
      myGraph % maxValence = maxValence
      myGraph % nColors    = 0

      ALLOCATE( myGraph % valence(1:nDOF), & 
                myGraph % color(1:nDOF), &
                myGraph % neighbors(1:maxValence,1:nDOF) )

      myGraph % valence    = 0
      myGraph % color      = 0
      myGraph % neighbors = 0

 END SUBROUTINE Build_POP_AdjacencyGraph
!
 SUBROUTINE Trash_POP_AdjacencyGraph( myGraph )
 ! This subroutine frees memory held by the allocatable attributes of the
 ! POP_AdjacencyGraph data structure.
 ! ======================================================================= !
   IMPLICIT NONE
   CLASS( POP_AdjacencyGraph ), INTENT(inout) :: myGraph

      DEALLOCATE( myGraph % valence, &
                  myGraph % color, &
                  myGraph % neighbors )

 END SUBROUTINE Trash_POP_AdjacencyGraph
!
!
!
 SUBROUTINE ConstructFromStencil_POP_AdjacencyGraph( myGraph, mesh, relStencil )
 ! Given a POP_Mesh and a Stencil, this routine constructs an adjacency
 ! graph using only the wet points in the mesh.

    IMPLICIT NONE
    CLASS( POP_AdjacencyGraph ), INTENT(inout) :: myGraph
    TYPE( POP_Mesh ), INTENT(in)               :: mesh
    TYPE( Stencil ), INTENT(in)                :: relStencil
    ! Local
    INTEGER :: i, j, k, e1, e2, m
    INTEGER :: this_i, this_j, this_k, true_i, true_j
    REAL(prec) :: t1, t2

       PRINT*, ' S/R ConstructFromStencil : Start! '

       CALL CPU_TIME( t1 )

       IF( .NOT.( ALLOCATED( myGraph % valence ) ) )THEN
          CALL myGraph % Build( mesh % nDOF, relStencil % nPoints ) 
       ENDIF


       DO e1 = 1, mesh % nDOF

          i = mesh % DOFtoIJK(1,e1)
          j = mesh % DOFtoIJK(2,e1)
          k = mesh % DOFtoIJK(3,e1)

          DO m = 1, relStencil % nPoints ! Loop over the stencil points
   
             ! Find the i,j,k indices for this point in the stencil around e1.
             this_i = i + relStencil % relativeNeighbors(1,m)
             this_j = j + relStencil % relativeNeighbors(2,m)
             this_k = k + relStencil % relativeNeighbors(3,m)
 
             ! The i and j components must be passed to a function for
             ! determining the actual i,j components -> This function call takes
             ! care of any periodicity/masking in the mesh.
             CALL GetTrueIJ( mesh % meshType, &
                             this_i, this_j, &
                             mesh % nX, mesh % nY, &
                             true_i, true_j )

             IF( true_i > 0 .AND. true_i <= mesh % nX .AND. &
                 true_j > 0 .AND. true_j <= mesh % nY .AND. &
                 this_k > 0 .AND. this_k <= mesh % nZ ) THEN

                ! Obtain the "degree of freedom" index for this neighbor 
                e2 = mesh % ijkToDOF( true_i, true_j, this_k )

                IF( e2 /= e1 .AND. e2 /= 0 )THEN
                   myGraph % valence(e1) = myGraph % valence(e1) + 1
                   myGraph % neighbors(myGraph % valence(e1), e1) = e2
                ENDIF

             ENDIF

          ENDDO

      ENDDO

      CALL CPU_TIME( t2 )
      PRINT *,' S/R ConstructFromStencil : Completed in ', t2-t1, ' seconds.'

 END SUBROUTINE ConstructFromStencil_POP_AdjacencyGraph
!
 SUBROUTINE GreedyColoring_POP_AdjacencyGraph( myGraph )
   ! Executes the Greedy coloring algorithm to fill in the "color" attribute
   ! given knowledge of the adjacency graph specified in the "neighbors"
   ! attribute.
   ! ====================================================================== !
   IMPLICIT NONE
   CLASS( POP_AdjacencyGraph ), INTENT(inout) :: myGraph
   ! Local
   INTEGER :: e1, e2, i, j
   INTEGER :: c2(1:myGraph % maxValence)
   LOGICAL :: neednewcolor, colorfound
   REAL(prec) :: t1, t2

      PRINT*, ' S/R GreedyColoring : Start! '
     
      CALL CPU_TIME( t1 )

      myGraph % ncolors = 1
      myGraph % color(1) = 1

      DO e1 = 2, myGraph % nDOF

         ! First build a list of all the colors of this element's neighbors
         c2 = 0
         DO i = 1, myGraph % valence(e1) 
            e2    = myGraph % neighbors(i,e1)
            c2(i) = myGraph % color(e2)
         ENDDO

         neednewcolor = .TRUE.

         DO i = 1, myGraph % nColors  ! only search through the colors already in use in the graph

           colorfound = .FALSE.
            DO j = 1, myGraph % valence(e1) ! Search through the neighbor-color list           
               IF( c2(j) == i )THEN         ! and determine if this color is already in use by any neighbors
                  colorfound = .TRUE.
               ENDIF
            ENDDO
            ! If we do not find color "i" then we can assign element "e1" the
            ! color "i"
            IF( .NOT. colorfound )THEN
               myGraph % color(e1) = i
               neednewcolor = .FALSE.
               EXIT
            ENDIF

         ENDDO

         ! In the event all of the possible colors are found in e1's neighbors,
         ! e1 must be colored with a new color      
         IF( neednewcolor )THEN
            myGraph % nColors   = myGraph % nColors + 1
            myGraph % color(e1) = myGraph % nColors
            PRINT*, 'nColors = ', myGraph % nColors
         ENDIF

      ENDDO
      
      CALL CPU_TIME( t2 )
      PRINT *,' S/R GreedyColoring : Coloring completed with ', myGraph % nColors, ' colors.'
      PRINT *,' S/R GreedyColoring : Coloring completed in ', t2-t1, ' seconds.'
      
 END SUBROUTINE GreedyColoring_POP_AdjacencyGraph
!
!
!
 SUBROUTINE ReadGraphFile_POP_AdjacencyGraph( myGraph, filename )
 !
   IMPLICIT NONE
   CLASS( POP_AdjacencyGraph ), INTENT(inout) :: myGraph
   CHARACTER(*), INTENT(in)                   :: filename
   ! Local
   INTEGER :: ndof, maxvalence, ncolors, i, j, fUnit

      OPEN( UNIT = NewUnit(fUnit), &
            FILE = TRIM(filename), &
            FORM = 'FORMATTED', &
            ACCESS = 'SEQUENTIAL', &
            STATUS = 'OLD', &
            ACTION = 'READ' )

      READ( fUnit, * ) ndof, maxvalence, ncolors

      CALL myGraph % Build( ndof, maxvalence )
      myGraph % nColors = ncolors
      
      DO i = 1, ndof

         READ( fUnit, * ) myGraph % valence(i)
         READ( fUnit, * ) myGraph % color(i)

         DO j = 1, myGraph % valence(i)
            READ( fUnit, * ) myGraph % neighbors(j,i)
         ENDDO
         
      ENDDO

      CLOSE( fUnit )

 END SUBROUTINE ReadGraphFile_POP_AdjacencyGraph
!
 SUBROUTINE WriteGraphFile_POP_AdjacencyGraph( myGraph, filename )
 !
   IMPLICIT NONE
   CLASS( POP_AdjacencyGraph ), INTENT(in) :: myGraph
   CHARACTER(*), INTENT(in)                :: filename
   ! Local
   INTEGER :: i, j, fUnit

      OPEN( UNIT = NewUnit(fUnit), &
            FILE = TRIM(filename), &
            FORM = 'FORMATTED', &
            ACCESS = 'SEQUENTIAL', &
            STATUS = 'REPLACE', &
            ACTION = 'WRITE' )

      WRITE( fUnit, * ) myGraph % nDOF, &
                        myGraph % maxValence, &
                        myGraph % nColors

      DO i = 1, myGraph % nDOF

         WRITE( fUnit, * ) myGraph % valence(i)
         WRITE( fUnit, * ) myGraph % color(i)

         DO j = 1, myGraph % valence(i)
            WRITE( fUnit, * ) myGraph % neighbors(j,i)
         ENDDO
         
      ENDDO

      CLOSE( fUnit )

 END SUBROUTINE WriteGraphFile_POP_AdjacencyGraph
!
 SUBROUTINE ReadGraphBinFile_POP_AdjacencyGraph( myGraph, filename )
 !
   IMPLICIT NONE
   CLASS( POP_AdjacencyGraph ), INTENT(inout) :: myGraph
   CHARACTER(*), INTENT(in)                   :: filename
   ! Local
   INTEGER :: ndof, maxvalence, ncolors, i, j, k, fUnit, datalength, reclength
   INTEGER :: nIntPerChunk, nChunks, adr1, adr2
   INTEGER, ALLOCATABLE :: localIOarray(:)

      OPEN( UNIT = NewUnit(fUnit), &
            FILE = TRIM(filename)//'.graph.hdr', &
            FORM = 'FORMATTED', &
            ACCESS = 'SEQUENTIAL', &
            STATUS = 'OLD', &
            ACTION = 'READ' )

      READ( fUnit, * ) ndof, maxvalence, ncolors

      CLOSE( UNIT = fUnit )

      CALL myGraph % Build( ndof, maxvalence )
      myGraph % nColors = ncolors
      
      datalength = ndof*(maxvalence+2)
      ! NOTE : in gfortran, internally the record length is stored as a 32-bit
      ! signed integer. Thus, in gfortran, it is not possible to have a record
      ! length be greater than 2^31 Bytes ( ~2.147 GB ). Our solution is to make
      ! the I/O array an integer multiple of a fixed "data chunk" size. For
      ! example, this chunk size might be 100 MB, which corresponds to
      ! 25-million 4-byte integers per record in the fil
      !
      ! Make sure that the datalength is an integer multiple of the
      ! "lgFileIOChunkSize". This parameter has units of "Bytes"
      IF( SIZEOF(i) == 4 )THEN

         ! Make sure that datalength is an integer multiple of nInt4PerChunk 
         nChunks    = (datalength/nInt4PerChunk)
         IF( nChunks*nInt4PerChunk < datalength )THEN
            nChunks = nChunks+1
         ENDIF
         dataLength   = nChunks*nInt4PerChunk
         recLength    = nInt4PerChunk*4
         nIntPerChunk = nInt4PerChunk

      ELSEIF( SIZEOF(i) == 8)THEN

         ! Make sure that datalength is an integer multiple of nInt4PerChunk 
         nChunks    = (datalength/nInt8PerChunk)
         IF( nChunks*nInt8PerChunk < datalength )THEN
            nChunks = nChunks+1
         ENDIF
         dataLength = nChunks*nInt8PerChunk
         recLength  = nInt8PerChunk*8
         nIntPerChunk = nInt8PerChunk

      ELSE
         PRINT*, 'Module POP_AdjacencyGraph_Class.f90 : S/R WriteGraphBinFile '
         PRINT*, '     Invalid INTEGER SIZE. STOPPING'
         STOP
      ENDIF
     
      ALLOCATE( localIOarray(1:datalength) )
      localIOarray = 0   
   
      OPEN( UNIT = NewUnit(fUnit), &
            FILE = TRIM(filename)//'.graph.bin', &
            FORM = 'UNFORMATTED', &
            ACCESS = 'DIRECT', &
            STATUS = 'OLD', &
            ACTION = 'READ', &
            RECL = reclength )

      PRINT*, 'S/R ReadGraphBinFile : Done!'
      PRINT*, '        '//TRIM(filename)//'.graph.bin,  file size (GB):',REAL(reclength)*REAL(nChunks)/10.0**9
      DO i = 1, nChunks
         adr1 = 1 + (i-1)*(nIntPerChunk)
         adr2 = adr1 + nIntPerChunk - 1
         READ( UNIT=fUnit, REC=i ) localIOarray(adr1:adr2)
      ENDDO
      CLOSE( fUnit )

      k = 0 
      DO i = 1, myGraph % ndof
         
         k = k+1
         myGraph % valence(i) = localIOArray(k)
         k = k+1
         myGraph % color(i)   = localIOArray(k)

         DO j = 1, myGraph % maxValence
            k = k+1
            myGraph % neighbors(j,i) = localIOArray(k)
         ENDDO
         
      ENDDO

      DEALLOCATE( localIOarray )

 END SUBROUTINE ReadGraphBinFile_POP_AdjacencyGraph
!
 SUBROUTINE WriteGraphBinFile_POP_AdjacencyGraph( myGraph, filename )
 !
   IMPLICIT NONE
   CLASS( POP_AdjacencyGraph ), INTENT(inout) :: myGraph
   CHARACTER(*), INTENT(in)                   :: filename
   ! Local
   INTEGER :: ndof, maxvalence, ncolors, i, j, fUnit, reclength
   INTEGER(KIND=8) :: datalength, k, adr1, adr2, nChunks
   INTEGER :: nIntPerChunk
   INTEGER, ALLOCATABLE :: localIOarray(:)

      OPEN( UNIT = NewUnit(fUnit), &
            FILE = TRIM(filename)//'.graph.hdr', &
            FORM = 'FORMATTED', &
            ACCESS = 'SEQUENTIAL', &
            STATUS = 'REPLACE', &
            ACTION = 'WRITE' )

      WRITE( fUnit, * ) myGraph  % ndof, myGraph % maxvalence, myGraph % ncolors

      CLOSE( UNIT = fUnit )

      
      datalength = myGraph % ndof*(myGraph % maxvalence+2)
      ! NOTE : in gfortran, internally the record length is stored as a 32-bit
      ! signed integer. Thus, in gfortran, it is not possible to have a record
      ! length be greater than 2^31 Bytes ( ~2.147 GB ). Our solution is to make
      ! the I/O array an integer multiple of a fixed "data chunk" size. For
      ! example, this chunk size might be 100 MB, which corresponds to
      ! 25-million 4-byte integers per record in the fil
      !
      ! Make sure that the datalength is an integer multiple of the
      ! "lgFileIOChunkSize". This parameter has units of "Bytes"
      IF( SIZEOF(i) == 4 )THEN

         ! Make sure that datalength is an integer multiple of nInt4PerChunk 
         nChunks    = (datalength/nInt4PerChunk)
         IF( nChunks*nInt4PerChunk < datalength )THEN
            nChunks = nChunks+1
         ENDIF
         dataLength   = nChunks*nInt4PerChunk
         recLength    = nInt4PerChunk*4
         nIntPerChunk = nInt4PerChunk

      ELSEIF( SIZEOF(i) == 8)THEN

         ! Make sure that datalength is an integer multiple of nInt8PerChunk 
         nChunks    = (datalength/nInt8PerChunk)
         IF( nChunks*nInt8PerChunk < datalength )THEN
            nChunks = nChunks+1
         ENDIF
         dataLength = nChunks*nInt8PerChunk
         recLength  = nInt8PerChunk*8
         nIntPerChunk = nInt8PerChunk

      ELSE
         PRINT*, 'Module POP_AdjacencyGraph_Class.f90 : S/R WriteGraphBinFile '
         PRINT*, '     Invalid INTEGER SIZE. STOPPING'
         STOP
      ENDIF
      PRINT*, 'READ THIS LINE : ', nChunks, myGraph % ndof, datalength
      ALLOCATE( localIOarray(1:datalength) )
      localIOarray = 0   
   
      k = 0 
      DO i = 1, myGraph % ndof
         
         k = k+1
         localIOArray(k) = myGraph % valence(i)
         k = k+1
         localIOArray(k) = myGraph % color(i)

         DO j = 1, myGraph % maxValence
            k = k+1
            localIOArray(k) = myGraph % neighbors(j,i)
         ENDDO
         
      ENDDO

      OPEN( UNIT = NewUnit(fUnit), &
            FILE = TRIM(filename)//'.graph.bin', &
            FORM = 'UNFORMATTED', &
            ACCESS = 'DIRECT', &
            STATUS = 'REPLACE', &
            ACTION = 'WRITE', &
            RECL = reclength )

      PRINT*, TRIM(filename)//'.graph.bin,  file size (GB):',REAL(reclength)*REAL(nChunks)/10.0**9
      DO i = 1, nChunks
         adr1 = 1 + (i-1)*(nIntPerChunk)
         adr2 = adr1 + nIntPerChunk - 1
         WRITE( UNIT=fUnit, REC=i ) localIOarray(adr1:adr2)
      ENDDO
      CLOSE( fUnit )

      DEALLOCATE( localIOarray )

 END SUBROUTINE WriteGraphBinFile_POP_AdjacencyGraph
END MODULE POP_AdjacencyGraph_Class
