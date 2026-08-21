SUBROUTINE PUT_SSTATS(stats, IPAR)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2007
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! write NetCDF output files -- summary statistics
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE multistats_types, only: SUMMARY                   ! summary statistics type
USE model_defn                                        ! model definition structures (includes filename)
USE meta_stats                                        ! metadata for summary statistics
USE model_numerix                                     ! model numerix parameters and arrays
USE sumextract_module                                 ! module to extract summary statistics
 USE handle_err_module, only: handle_err              ! handle NetCDF errors
IMPLICIT NONE
! input
type(SUMMARY), intent(in)              :: stats       ! structures that depend on nState/nPar
INTEGER(I4B),    INTENT(IN)            :: IPAR        ! parameter set index
! internal
INTEGER(I4B)                           :: IERR,NCID   ! error code; NetCDF ID
INTEGER(I4B), DIMENSION(1)             :: INDX        ! indices for parameter write
INTEGER(I4B)                           :: IVAR        ! loop through parameters
REAL(WP)                               :: XPAR        ! desired parameter (working precision)
REAL(MSP)                              :: APAR        ! desired parameter (MSP is SP)
INTEGER(I4B)                           :: IVAR_ID     ! variable ID
include 'netcdf.inc'                                  ! use netCDF libraries
! ---------------------------------------------------------------------------------------
! open NetCDF parameter file
IERR = NF_OPEN(TRIM(FNAME_NETCDF_PARA),NF_WRITE,NCID); CALL HANDLE_ERR(IERR)

 ! define indices for model output
 INDX = (/IPAR/)

 ! loop through summary statistics
 DO IVAR=1,NSUMVAR
  XPAR = SUMEXTRACT( stats, XNAME(IVAR) ); APAR=XPAR                         ! get parameter ivar
  IERR = NF_INQ_VARID(NCID,TRIM(XNAME(IVAR)),IVAR_ID); CALL HANDLE_ERR(IERR) ! get variable ID
  IERR = NF_PUT_VAR1_REAL(NCID,IVAR_ID,INDX,APAR); CALL HANDLE_ERR(IERR)     ! write data
 END DO  ! (ivar)

IERR = NF_CLOSE(NCID)
! ---------------------------------------------------------------------------------------
END SUBROUTINE PUT_SSTATS
