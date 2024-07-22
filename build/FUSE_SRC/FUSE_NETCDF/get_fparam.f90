module get_fparam_module
USE nrtype
USE netcdf
implicit none
private
public::GET_PRE_PARAM,GET_SCE_PARAM

contains

SUBROUTINE GET_PRE_PARAM(NETCDF_FILE,ISET,IMOD,MPAR,XPAR)

! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Nans Addor, 2017 - Based on GET_SCE_PARAM
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Load a pre-defined set of parameter values from a NetCDF file
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE fuse_fileManager, only : OUTPUT_PATH              ! define output path
USE multiparam, ONLY: LPARAM, NUMPAR                  ! parameter names
IMPLICIT NONE
! input
CHARACTER(LEN=*), INTENT(IN)           :: NETCDF_FILE ! NetCDF file name
INTEGER(I4B), INTENT(IN)               :: ISET        ! indice of parameter set to extract
INTEGER(I4B), INTENT(IN)               :: IMOD        ! model index
INTEGER(I4B), INTENT(IN)               :: MPAR        ! number of model parameters
! internal
INTEGER(I4B), DIMENSION(1)             :: INDX        ! indices for parameter extraction
LOGICAL(LGT)                           :: LEXIST      ! .TRUE. if NetCDF file exists
INTEGER(I4B)                           :: IERR        ! error code
INTEGER(I4B)                           :: NCID        ! NetCDF file ID
INTEGER(I4B)                           :: IDIMID      ! NetCDF dimension ID
INTEGER(I4B)                           :: IVARID      ! NetCDF variable ID
INTEGER(I4B)                           :: IPAR        ! loop through model parameters
INTEGER(I4B)                           :: NPAR        ! number of parameter sets in output file
REAL(DP)                               :: APAR        ! parameter value (single precision)
! output
REAL(SP), DIMENSION(MPAR), INTENT(OUT) :: XPAR        ! parameter value (whatever precision SP is)
include 'netcdf.inc'                                  ! use netCDF libraries
! ---------------------------------------------------------------------------------------
! check that the file exists
INQUIRE(FILE=TRIM(NETCDF_FILE),EXIST=LEXIST)
IF (.NOT.LEXIST) THEN
 print *, 'The NetCDF file containing the predefined parameter set does not exist:'
 print *, TRIM(NETCDF_FILE)
 STOP
ENDIF

print *, 'Opening parameter file:', TRIM(NETCDF_FILE)

! open file
IERR = NF_OPEN(TRIM(NETCDF_FILE),NF_NOWRITE,NCID); CALL HANDLE_ERR(IERR)

! get number of parameter sets
IERR = NF_INQ_DIMID(NCID,'par',IDIMID); CALL HANDLE_ERR(IERR)
IERR = NF_INQ_DIMLEN(NCID,IDIMID,NPAR); CALL HANDLE_ERR(IERR)

print *, 'NPAR - total number of parameter sets in parameter file:', NPAR

IF (ISET.GT.NPAR) THEN
 print *, 'Impossible to extract parameter set', ISET, 'since there are only', NPAR, 'parameter sets'
 STOP
ENDIF

print *, 'Extracting parameter set', ISET

! loop through parameters
DO IPAR=1,NUMPAR

  ! get parameter id
  IERR = NF_INQ_VARID(NCID,TRIM(LPARAM(IPAR)%PARNAME),IVARID); CALL HANDLE_ERR(IERR)

  ! get parameter value for the selected parameter set
  INDX = (/ISET/)
  IERR = NF_GET_VAR1_DOUBLE(NCID,IVARID,INDX,APAR); CALL HANDLE_ERR(IERR)

  ! put parameter value in the output vector
  XPAR(IPAR) = APAR

  print *, 'PARAM VALUES:',LPARAM(IPAR)%PARNAME, '->', APAR

END DO

PRINT *, 'Predefined parameter set loaded into XPAR!'

! close NetCDF file
IERR = NF_CLOSE(NCID)
! ---------------------------------------------------------------------------------------
END SUBROUTINE GET_PRE_PARAM

SUBROUTINE GET_SCE_PARAM(NETCDF_FILE,IMOD,MPAR,XPAR)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2009
! Modified by Nans Addor - previously the parameter set extracted was the last one,
! now it is the one associated with the lowest RMSE
! Modified by Cyril Thébault to add KGE metric, 2024
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Read parameters in LPARAM from the parameter set with the highest KGE in the specified NetCDF file
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE fuse_fileManager, only : OUTPUT_PATH              ! define output path
USE multiparam, ONLY: LPARAM, NUMPAR                  ! parameter names
IMPLICIT NONE
! input
CHARACTER(LEN=*), INTENT(IN)           :: NETCDF_FILE ! NetCDF file name
INTEGER(I4B), INTENT(IN)               :: IMOD        ! model index
INTEGER(I4B), INTENT(IN)               :: MPAR        ! number of model parameters
! internal
INTEGER(I4B), DIMENSION(1)             :: INDX        ! indices for parameter extraction
LOGICAL(LGT)                           :: LEXIST      ! .TRUE. if NetCDF file exists
INTEGER(I4B)                           :: IERR        ! error code
INTEGER(I4B)                           :: NCID        ! NetCDF file ID
INTEGER(I4B)                           :: IDIMID      ! NetCDF dimension ID
INTEGER(I4B)                           :: IVARID      ! NetCDF variable ID
INTEGER(I4B)                           :: I_KGE       ! NetCDF KGE ID
INTEGER(I4B)                           :: I_OPT_PARA  ! index of the optimum parameter set (e.g. lowest RSME) - MUST BE DIMENSIONED
REAL(DP), DIMENSION(:),ALLOCATABLE     :: KGE         ! KGE for each parameter set
INTEGER(I4B), DIMENSION(1)             :: ARRAY_SIZE  ! Dimension of the KGE 
REAL(DP)                               :: HIGHEST_KGE ! Highest KGE
INTEGER(I4B)                           :: IPAR        ! loop through model parameters
INTEGER(I4B)                           :: NPAR        ! number of parameter sets in output file
REAL(DP)                               :: APAR        ! parameter value (single precision)
! output
REAL(SP), DIMENSION(MPAR), INTENT(OUT) :: XPAR        ! parameter value (whatever precision SP is)
include 'netcdf.inc'                                  ! use netCDF libraries
! ---------------------------------------------------------------------------------------
! check that the file exists
INQUIRE(FILE=TRIM(NETCDF_FILE),EXIST=LEXIST)
IF (.NOT.LEXIST) THEN
 print *, 'The NetCDF file containing the SCE parameter sets does not exist:'
 print *, TRIM(NETCDF_FILE)
 STOP
ELSE

  print *, 'Loading SCE parameters from'
  print *, TRIM(NETCDF_FILE)

ENDIF

! open file
IERR = NF_OPEN(TRIM(NETCDF_FILE),NF_NOWRITE,NCID); CALL HANDLE_ERR(IERR)

 ! get number of parameter sets
 IERR = NF_INQ_DIMID(NCID,'par',IDIMID); CALL HANDLE_ERR(IERR)
 IERR = NF_INQ_DIMLEN(NCID,IDIMID,NPAR); CALL HANDLE_ERR(IERR)

 ! extract KGE for each parameter set
 print *, 'Length of the par dimension (the number of parameter sets produced by SCE is higher)', NPAR

 ALLOCATE(KGE(NPAR),STAT=IERR); IF(IERR.NE.0) STOP ' problem allocating space for KGE '

 IERR = NF_INQ_VARID(NCID,'kge',I_KGE); CALL HANDLE_ERR(IERR)
 IERR = NF_GET_VAR_DOUBLE(NCID,I_KGE,KGE); CALL HANDLE_ERR(IERR)
 
 ! Find the best value lower than opt value
 HIGHEST_KGE = -9999
 I_OPT_PARA = -9999
 
 DO IPAR = 1, NPAR
   IF (KGE(IPAR) .LE. 9.0E+35 .AND. KGE(IPAR) .GT. HIGHEST_KGE) THEN
     HIGHEST_KGE = KGE(IPAR)
     I_OPT_PARA = IPAR
   END IF
 END DO

 print *, 'Index of parameter set with highest KGE =',I_OPT_PARA
 print *, 'Highest KGE =',HIGHEST_KGE

 DEALLOCATE(KGE,STAT=IERR); IF (IERR.NE.0) STOP ' problem deallocating ATIME/TDATA '

 PRINT *, 'Reading from NetCDF file parameter values for best parameter set:'

 ! loop through parameters
 DO IPAR=1,NUMPAR

  ! get parameter id
  IERR = NF_INQ_VARID(NCID,TRIM(LPARAM(IPAR)%PARNAME),IVARID); CALL HANDLE_ERR(IERR)

  ! get parameter value for the optimal parameter set
  INDX = (/I_OPT_PARA/)
  IERR = NF_GET_VAR1_DOUBLE(NCID,IVARID,INDX,APAR); CALL HANDLE_ERR(IERR)

  ! put parameter value in the output vector
  XPAR(IPAR) = APAR

  print *, 'PARAM VALUES:',LPARAM(IPAR)%PARNAME, '->', APAR

 END DO

 PRINT *, 'Best parameter set loaded into XPAR!'

! close NetCDF file
IERR = NF_CLOSE(NCID)
! ---------------------------------------------------------------------------------------
END SUBROUTINE GET_SCE_PARAM

end module get_fparam_module
