! CRSMatrix_Class.f90
! 
! Copyright 2016 Joseph Schoonover, Los Alamos National Laboratory (jschoonover@lanl.gov) 
! All rights reserved. 
! 
! CRSMatrix_Class.f90 is part of the Fast Equilibration of Ocean Tracers Software (FEOTS). 
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !
 
MODULE CRSMatrix_Class
!
! This module contains a data-structure and associated routines for working with sparse matrices.
! We use the "Compressed Row Storage" model for handling sparse matrix operations. In this model,
! the elements of a matrix are stored in a single dimension array of real values. The indices
! for the start and end of each row are stored in a 2-D integer array and the columns of each
! nonzero entry are stored. As an example, consider the tri-diagonal matrix
!
! ______________________________
! Column : 1  2  3  4  5  6  7  | row |
! ______________________________|     |
!         -2  1  0  0  0  0  0  |  1  |
!          1 -2  1  0  0  0  0  |  2  |
!          0  1 -2  1  0  0  0  |  3  |
!          0  0  1 -2  1  0  0  |  4  |
!          0  0  0  1 -2  1  0  |  5  |
!          0  0  0  0  1 -2  1  |  6  |
!          0  0  0  0  0  1 -2  |  7  |
!
!  There are 7 + 6 + 6 = 19 elements in this 7x7 tri-diagonal matrix.
!  In the CRS structure, the matrix data would be stored as
!
!  A(1) = -2, A(2) = 1,
!  A(3) = 1,  A(4) = -2, A(5) = 1,
!  A(6) = 1,  A(7) = -2, A(8) = 1,
!   ........
!  A(18) = 1, A(19) = -2
!
!  The first row contains elements 1 and 2, second row contains elements (3,4,5), etc. We would 
!  store the start and end of the rows as follows :
!
!  rowBounds(1,1) = 1, rowBounds(1,2) = 2
!  rowBounds(2,1) = 3, rowBounds(2,2) = 5
!  rowBounds(3,1) = 6, rowBounds(3,2) = 8
!   ........
!  rowBounds(7,1) = 18, rowBounds(7,2) = 19
!
!  * Notice that the first index of "rowBounds" varies from 1 to the number of rows.
!
!  Finally, the column needs to be stored for each non-zero element.
!
!  col(1) = 1, col(2) = 2
!  col(3) = 1, col(4) = 2, col(5) = 3
!  col(6) = 2, col(7) = 3, col(8) = 4
!   .......
!  col(18) = 6, col(19) = 7
!
! ================================================================================================ !

 USE ConstantsDictionary
 USE CommonRoutines

 IMPLICIT NONE


   TYPE CRSMatrix  
      INTEGER                         :: nRows, nCols, nElems
      REAL(prec), ALLOCATABLE         :: A(:)
      INTEGER, ALLOCATABLE            :: rowBounds(:,:), col(:)

      CONTAINS

      PROCEDURE :: Build => Build_CRSMatrix
      PROCEDURE :: Reset => Reset_CRSMatrix
      PROCEDURE :: Trash => Trash_CRSMatrix

      PROCEDURE :: GetDataForRowAndCol => GetDataForRowAndCol_CRSMatrix
      PROCEDURE :: MatVecMul           => MatVecMul_CRSMatrix
      PROCEDURE :: DenseToSparse       => DenseToSparse_CRSMatrix
      PROCEDURE :: SubSample           => SubSample_CRSMatrix

!      PROCEDURE :: ReadHeader  => ReadHeader_CRSMatrix
      PROCEDURE :: WriteHeader             => WriteHeader_CRSMatrix
      PROCEDURE :: ReadSparseConnectivity  => ReadSparseConnectivity_CRSMatrix
      PROCEDURE :: WriteSparseConnectivity => WriteSparseConnectivity_CRSMatrix
      PROCEDURE :: ReadMatrixData          => ReadMatrixData_CRSMatrix
      PROCEDURE :: WriteMatrixData         => WriteMatrixData_CRSMatrix

   END TYPE CRSMatrix

   

 CONTAINS
!
!
!==================================================================================================!
!------------------------------- Manual Constructors/Destructors ----------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE Build_CRSMatrix( myMatrix, nRows, nCols, nElems )
 ! S/R Build
 !
 !  
 ! =============================================================================================== !
  IMPLICIT NONE
  CLASS( CRSMatrix ), INTENT(out) :: myMatrix
  INTEGER, INTENT(in)             :: nRows, nCols, nElems

     myMatrix % nRows = nRows
     myMatrix % nCols = nCols
     myMatrix % nElems = nElems

     ALLOCATE( myMatrix % A(1:nElems), &
               myMatrix % rowBounds(1:2,1:nRows), &
               myMatrix % col(1:nElems) )

     myMatrix % A         = 0.0_prec
     myMatrix % rowBounds = 0
     myMatrix % col       = 0

 END SUBROUTINE Build_CRSMatrix
!
 SUBROUTINE Reset_CRSMatrix( myMatrix )
 ! S/R Reset
 !
 !  
 ! =============================================================================================== !
  IMPLICIT NONE
  CLASS( CRSMatrix ), INTENT(inout) :: myMatrix

     myMatrix % A         = 0.0_prec
     myMatrix % rowBounds = 0
     myMatrix % col       = 0

 END SUBROUTINE Reset_CRSMatrix
!
 SUBROUTINE Trash_CRSMatrix( myMatrix )
 ! S/R Trash
 !
 !  
 ! =============================================================================================== !
  IMPLICIT NONE
  CLASS( CRSMatrix ), INTENT(inout) :: myMatrix

     DEALLOCATE( myMatrix % rowBounds, myMatrix % col, myMatrix % A )
  
 END SUBROUTINE Trash_CRSMatrix
!
!
!==================================================================================================!
!--------------------------------- Type Specific Routines  ----------------------------------------!
!==================================================================================================!
!
!
 FUNCTION GetDataForRowAndCol_CRSMatrix( myMatrix, i, j ) RESULT( outData )
 ! S/R GetDataForRowAndCol
 !
 ! =============================================================================================== !
   CLASS( CRSMatrix ) :: myMatrix
   REAL(prec)         :: outData
   INTEGER            :: i, j
   ! LOCAL
   INTEGER :: rstart, rEnd, k, col
   INTEGER :: localCol(1:myMatrix % nElems)
   
      outData = 0.0_prec
      
      localCol = myMatrix % col

      rstart = myMatrix % rowBounds(1,i)
      rEnd = myMatrix % rowBounds(2,i)

      k = rstart
      col = 0
      DO WHILE( col /= j .AND. k <= rEnd )
         col = localCol(k)
         k = k + 1
      ENDDO

      IF( k <= rEnd )THEN
         outData = myMatrix % A(k)
      ENDIF

 END FUNCTION GetDataForRowAndCol_CRSMatrix
!
!
!
 FUNCTION MatVecMul_CRSMatrix( myMatrix, x ) RESULT( Ax )
 ! FUNCTION MatVecMul
 !
 !    Performs matrix-vector multiply with the sparse matrix and a vector (stored as 1-D array).
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( CRSMatrix ) :: myMatrix
   REAL(prec)         :: x(1:myMatrix % nCols)
   REAL(prec)         :: Ax(1:myMatrix % nRows)
   ! Local
   INTEGER    :: iel, row
   REAL(prec) :: rowsum

      ! For each row of the matrix, a vector-product must be computed; recall that a matrix-vector
      ! product can be written :
      !                          (Ax)_i = sum_{j}( A_{i,j}x_{j} )
      ! This is equivalent to
      !                          (Ax)_i = sum_{j,A_{i,j} /= 0.0}( A_{i,j}x_{j} )
      ! where the sum is taken over non-zero elements. To achieve this, we cycle over the rows and
      ! find the non-zero columns (j) and store those input vector values. 

      
      ! Then, we cycle over the rows and compute the sums over each row.
      DO row = 1, myMatrix % nRows
         
         rowsum = 0.0_prec
         DO iel = myMatrix % rowBounds(1,row), myMatrix % rowBounds(2,row)
            rowSum = rowSum + myMatrix % A(iel)*x(myMatrix % col(iel))
         ENDDO
         Ax(row) = rowSum

      ENDDO
      
 END FUNCTION MatVecMul_CRSMatrix
!
!
!
 SUBROUTINE DenseToSparse_CRSMatrix( myMatrix, dmatrix, nRows, nCols )
 ! S/R DenseToSparse
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( CRSMatrix ), INTENT(out) :: myMatrix
   INTEGER, INTENT(in)             :: nRows, nCols
   REAL(prec), INTENT(in)          :: dMatrix(1:nRows,1:nCols)
   ! Local
   INTEGER    :: row, col, nEl, r1, r2
   REAL(prec) :: aij


      ! First, count the number of non-zero entries.
      nEl = 0
      DO col = 1, nCols
         DO row = 1, nRows
         
            aij = dMatrix(row,col)

            IF( .NOT.( AlmostEqual(aij, ZERO) ) )THEN
               nEl = nEl+1
            ENDIF

         ENDDO
      ENDDO
      PRINT*, nEl
      CALL myMatrix % Build( nRows, nCols, nEl )

      nEl = 0
      r1 = 0
      DO row = 1, nRows
         DO col = 1, nCols
         
            aij = dMatrix(row,col)

            IF( .NOT.( AlmostEqual(aij, ZERO) ) )THEN
               nEl = nEl+1
               myMatrix % A(nEl)   = aij
               myMatrix % col(nEl) = col
               ! Check to see if we've changed rows
               r2 = row
               IF( r2 /= r1 )THEN
                  myMatrix % rowBounds(1,row) = nEl
               ENDIF
               r1 = r2 
            ENDIF

         ENDDO
      ENDDO

      DO row = 2, nRows
         myMatrix % rowBounds(2,row-1) = myMatrix % rowBounds(1,row)-1
      ENDDO
      myMatrix % rowBounds(2,nRows) = nEl
      

 END SUBROUTINE DenseToSparse_CRSMatrix
!
!
 SUBROUTINE SubSample_CRSMatrix( inputMatrix, subMatrix, dofMap, inverseDOFMap, nDOF, nMaxElem )
! This subroutine constructs a sparse matrix that contains only rows and
! columns specified in the "dofMap" - it "subsamples" the matrix.
   IMPLICIT NONE
   CLASS( CRSMatrix ), INTENT(in)   :: inputMatrix
   TYPE( CRSMatrix ), INTENT(inout) :: subMatrix
   INTEGER, INTENT(in)              :: nDOF, nMaxElem
   INTEGER, INTENT(in)              :: dofMap(1:nDOF)
   INTEGER, INTENT(in)              :: inverseDOFMap(1:inputMatrix % nRows)
   ! Local
   INTEGER :: regionMask(1:inputMatrix % nRows)
   INTEGER :: i, j, iEl, row, col, nInRow, e1, e2, jlocal
   REAL(prec) :: rowdata(1:20)

      IF( .NOT. ALLOCATED( subMatrix % A ) )THEN
         ! Allocate space for the submatrix
         ALLOCATE( subMatrix % A(1:nMaxElem), &
                   subMatrix % rowBounds(1:2,1:nDOF), &
                   subMatrix % col(1:nMaxElem) )
         subMatrix % nRows  = nDOF
         subMatrix % nCols  = nDOF
         subMatrix % nElems = nMaxElem
      ENDIF
      ! To subsample, the sparse matrix, we first need to determine the number
      ! of non-zero elements contained within the dofMap
      !  The "dofMap" tells us which rows to keep. We can cycle over each row of
      !  the sparse matrix and check a regional mask to determine if the column
      !  value should be kept.
      ! The region mask can be constructed using the dofMap. First all of the
      ! entries of the region mask are set to zero.
      regionMask = 0
      ! Then, we cycle through the dofMap to "turn on" points only set in the
      ! dofMap.
      DO i = 1, nDOF
         regionMask( dofMap(i) ) = 1
      ENDDO

      rowdata = 0.0_prec
      ! Now we fill in the submatrix data
      iEl = 0
      subMatrix % rowBounds(1,1) = 1
      DO i = 1, nDOF
         e1     = inputMatrix % rowBounds(1,dofMap(i))
         e2     = inputMatrix % rowBounds(2,dofMap(i))
         nInRow = e2-e1 + 1

         rowdata(1:nInRow) = inputMatrix % A(e1:e2)

         jlocal = 0
         DO j = e1, e2
            jlocal = jlocal + 1
            IF( regionMask(inputMatrix % col(j)) == 1 )THEN
               iEl = iEl + 1
               subMatrix % A(iEl)   = rowdata( jlocal )
               subMatrix % col(iEl) = inverseDOFMap( inputMatrix % col(j) ) ! returns the "regional" degree of freedom for the global degree of freedom 
            ENDIF
         ENDDO

         subMatrix % rowBounds( 2, inverseDOFMap( dofMap(i) ) ) = iEl
         IF( inverseDOFMap( dofMap(i) ) + 1 <= nDOF )THEN
            subMatrix % rowBounds( 1, inverseDOFMap( dofMap(i) )+1 ) = iEl+1
         ENDIF

      ENDDO

 END SUBROUTINE SubSample_CRSMatrix
!
!
!==================================================================================================!
!-------------------------------------- File I/O Routines  ----------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE ReadHeader( fileBase, nRow, nCol, nElems )
 ! S/R ReadHeader
 !
 !   
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CHARACTER(*), INTENT(in)          :: fileBase
   INTEGER, INTENT(out)              :: nRow, nCol, nElems
   ! Local
   INTEGER :: fUnit

      OPEN( UNIT = NewUnit(fUnit), & 
            FILE = TRIM(fileBase)//'.head', &
            FORM = 'FORMATTED', &
            ACCESS = 'SEQUENTIAL', &
            STATUS = 'OLD', &
            ACTION = 'READ' )

      READ( fUnit, * ) nrow, ncol, nElems

      CLOSE( fUnit )

 END SUBROUTINE ReadHeader
!
 SUBROUTINE WriteHeader_CRSMatrix( myMatrix, fileBase )
 ! S/R WriteHeader
 !
 !   
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( CRSMatrix ), INTENT(inout) :: myMatrix
   CHARACTER(*), INTENT(in)          :: fileBase
   ! Local
   INTEGER :: fUnit

      OPEN( UNIT = NewUnit(fUnit), & 
            FILE = TRIM(fileBase)//'.head', &
            FORM = 'FORMATTED', &
            ACCESS = 'SEQUENTIAL', &
            STATUS = 'REPLACE', &
            ACTION = 'WRITE' )

      WRITE( fUnit, * ) myMatrix % nRows, myMatrix % nCols, myMatrix % nElems

      CLOSE( fUnit )

 END SUBROUTINE WriteHeader_CRSMatrix
!
 SUBROUTINE ReadSparseConnectivity_CRSMatrix( myMatrix, fileBase )
 ! S/R ReadSparseConnectivity
 !
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( CRSMatrix ), INTENT(inout) :: myMatrix
   CHARACTER(*), INTENT(in)          :: fileBase
   ! Local
   INTEGER :: i
   INTEGER :: recID
   INTEGER :: fUnit


      OPEN( UNIT = NewUnit(fUnit), & 
            FILE = TRIM(fileBase)//'.conn', &
            FORM = 'UNFORMATTED', &
            ACCESS = 'DIRECT', &
            STATUS = 'OLD', &
            ACTION = 'READ', &
            RECL = 4 )

      recID  = 1

      DO i = 1, myMatrix % nRows
         READ( fUnit, REC = recID ) myMatrix %  rowBounds(1,i)
         recID = recID + 1
         READ( fUnit, REC = recID ) myMatrix %  rowBounds(2,i)
         recID = recID + 1
      ENDDO
 
      DO i = 1, myMatrix % nElems
         READ( fUnit, REC=recID ) myMatrix % col(i)
         recID = recID + 1
      ENDDO

      CLOSE( fUnit )

      !DO i = 1, myMatrix % nRows-1
      !   myMatrix % rowBounds(i,2) = myMatrix % rowBounds(i+1,1)-1
      !ENDDO


 END SUBROUTINE ReadSparseConnectivity_CRSMatrix
!
 SUBROUTINE WriteSparseConnectivity_CRSMatrix( myMatrix, fileBase )
 ! S/R WriteSparseConnectivity
 !
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( CRSMatrix ), INTENT(inout) :: myMatrix
   CHARACTER(*), INTENT(in)          :: fileBase
   ! Local
   INTEGER :: i
   INTEGER :: recID
   INTEGER :: fUnit


      OPEN( UNIT = NewUnit(fUnit), & 
            FILE = TRIM(fileBase)//'.conn', &
            FORM = 'UNFORMATTED', &
            ACCESS = 'DIRECT', &
            STATUS = 'REPLACE', &
            ACTION = 'WRITE', &
            RECL = 4 )

      recID  = 1

      DO i = 1, myMatrix % nRows
         WRITE( fUnit, REC = recID ) myMatrix %  rowBounds(1,i)
         recID = recID + 1
         WRITE( fUnit, REC = recID ) myMatrix %  rowBounds(2,i)
         recID = recID + 1
      ENDDO
 
      DO i = 1, myMatrix % nElems
         WRITE( fUnit, REC=recID ) myMatrix % col(i)
         recID = recID + 1
      ENDDO

      CLOSE( fUnit )


 END SUBROUTINE WriteSparseConnectivity_CRSMatrix
!
 SUBROUTINE ReadMatrixData_CRSMatrix( myMatrix, fileBase )
 ! S/R ReadMatrixData
 !
 !    This subroutine reads the "fileBase.data" file that contains the matrix element data. As a 
 !    prerequisite to calling this routine, the "% nElems" attribute needs to be filled in and 
 !    memory for the "% A" attribute should be allocated.
 !    
 !    Usage :
 !       CALL myMatrix % ReadMatrixData( fileBase )
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( CRSMatrix ), INTENT(inout) :: myMatrix
   CHARACTER(*), INTENT(in)          :: fileBase
   ! Local
   INTEGER :: ioerr, nChunks, recLength, nRealPerChunk
   INTEGER :: fUnit, dataLength, s1, s2, i
   INTEGER :: row, iEl
   REAL(prec), ALLOCATABLE :: localIOarray(:)

       dataLength = myMatrix % nElems
       IF( prec==SP )THEN

         ! Make sure that datalength is an integer multiple of nInt4PerChunk 
         nChunks    = (datalength/nReal4PerChunk)
         IF( nChunks*nReal4PerChunk < datalength )THEN
            nChunks = nChunks+1
         ENDIF
         dataLength   = nChunks*nReal4PerChunk
         recLength    = nReal4PerChunk*prec
         nRealPerChunk = nReal4PerChunk

      ELSEIF( prec == DP )THEN

         ! Make sure that datalength is an integer multiple of nInt4PerChunk 
         nChunks    = (datalength/nReal8PerChunk)
         IF( nChunks*nReal8PerChunk < datalength )THEN
            nChunks = nChunks+1
         ENDIF
         dataLength = nChunks*nReal8PerChunk
         recLength  = nReal8PerChunk*prec
         nRealPerChunk = nReal8PerChunk

      ELSE
         STOP
      ENDIF

      ALLOCATE( localIOarray(1:datalength) )
      localIOarray = 0.0_prec
      !localIOarray(1:myMatrix % nElems) = myMatrix % A(1:myMatrix % nElems) 

      OPEN( UNIT = NewUnit(fUnit), & 
            FILE = TRIM(fileBase)//'.data', &
            FORM = 'UNFORMATTED', &
            ACCESS = 'DIRECT', &
            STATUS = 'OLD', &
            ACTION = 'READ', &
            RECL = recLength )
      
      DO i = 1, nChunks
         s1 = 1  + (nRealPerChunk)*(i-1)
         s2 = s1 + nRealPerChunk - 1
         READ( fUnit, REC=i, IOSTAT = ioerr ) localIOarray(s1:s2)
      ENDDO
      CLOSE( fUnit )

 
      myMatrix % A = 0.0_prec
      myMatrix % A(1:myMatrix % nElems) = localIOarray(1:myMatrix % nElems)
 
!      DO row = 1, myMatrix % nRows
!         
!         DO iel = myMatrix % rowBounds(1,row), myMatrix % rowBounds(2,row)
!            IF( ABS( myMatrix % A(iel) ) > 4.0_prec*10.0_prec**(-4) )THEN
!              PRINT*, myMatrix % A(iel), row, myMatrix % col(iel)
!            ENDIF
!         ENDDO
!
!      ENDDO
      
     
      DEALLOCATE( localIOarray )

 END SUBROUTINE ReadMatrixData_CRSMatrix
!
 SUBROUTINE WriteMatrixData_CRSMatrix( myMatrix, fileBase )
 ! S/R WriteMatrixData
 !
 !    Usage :
 !       CALL myMatrix % WriteMatrixData( fileBase )
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( CRSMatrix ), INTENT(in) :: myMatrix
   CHARACTER(*), INTENT(in)       :: fileBase
   ! Local
   INTEGER :: ioerr, nChunks, recLength, nRealPerChunk
   INTEGER :: fUnit, dataLength, s1, s2, i
   REAL(prec), ALLOCATABLE :: localIOarray(:)

       dataLength = myMatrix % nElems
       IF( prec==SP )THEN

         ! Make sure that datalength is an integer multiple of nInt4PerChunk 
         nChunks    = (datalength/nReal4PerChunk)
         IF( nChunks*nReal4PerChunk < datalength )THEN
            nChunks = nChunks+1
         ENDIF
         dataLength   = nChunks*nReal4PerChunk
         recLength    = nReal4PerChunk*prec
         nRealPerChunk = nReal4PerChunk

      ELSEIF( prec == DP )THEN

         ! Make sure that datalength is an integer multiple of nInt4PerChunk 
         nChunks    = (datalength/nReal8PerChunk)
         IF( nChunks*nReal8PerChunk < datalength )THEN
            nChunks = nChunks+1
         ENDIF
         dataLength = nChunks*nReal8PerChunk
         recLength  = nReal8PerChunk*prec
         nRealPerChunk = nReal8PerChunk

      ELSE
         STOP
      ENDIF

      ALLOCATE( localIOarray(1:datalength) )
      localIOarray = 0.0_prec
      localIOarray(1:myMatrix % nElems) = myMatrix % A(1:myMatrix % nElems) 

      OPEN( UNIT = NewUnit(fUnit), & 
            FILE = TRIM(fileBase)//'.data', &
            FORM = 'UNFORMATTED', &
            ACCESS = 'DIRECT', &
            STATUS = 'REPLACE', &
            ACTION = 'WRITE', &
            RECL = recLength )
      
      DO i = 1, nChunks
         s1 = 1  + (nRealPerChunk)*(i-1)
         s2 = s1 + nRealPerChunk - 1
         WRITE( fUnit, REC=i, IOSTAT = ioerr ) localIOarray(s1:s2)
      ENDDO
      CLOSE( fUnit )

      DEALLOCATE( localIOarray )

 END SUBROUTINE WriteMatrixData_CRSMatrix
!
!
!==================================================================================================!
!------------------------------------ File I/O Routines  ------------------------------------------!
!==================================================================================================!
!
!
END MODULE CRSMatrix_Class
