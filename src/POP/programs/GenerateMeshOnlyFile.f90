PROGRAM GenerateMeshOnlyFile

! src/common/
USE CommonRoutines
! src/POP/
USE POP_Params_Class
USE POP_Mesh_Class

IMPLICIT NONE

   TYPE( POP_Params ) :: params
   TYPE( POP_Mesh )   :: mesh
   INTEGER(KIND=8), ALLOCATABLE :: dofToIJK_check(:,:), ijkToDOF_check(:,:,:)
   INTEGER :: i, j, k, fUnit, diffcount, thisdiff
   CHARACTER(400) :: ncfile


      CALL params % Build( )

      ! Reads in the first file from the IRF File list
      OPEN( UNIT=NewUnit(fUnit),&
            FILE=TRIM(params % IRFListFile), &
            FORM='FORMATTED',&
            ACCESS='SEQUENTIAL',&
            ACTION='READ',&
            STATUS='OLD' )
      READ( fUnit, '(A400)' ) ncFile
      CLOSE( fUnit )

      ! This call loads the mesh from the netcdf file and builds the
      ! ijkToDOF and dofToIJK mappings by using the "kmt" field.
      CALL mesh % Load( TRIM(ncFile)  )

      ! Write a netcdf file containing only the mesh
      CALL mesh % WriteNetCDF( TRIM(params % meshFile) )


      CALL mesh % Trash( )

END PROGRAM GenerateMeshOnlyFile

