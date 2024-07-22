FUNCTION FUNCTN(NOPT,A)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2009
! Modified by Cyril Thébault to add KGE metric, 2024
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Wrapper for SCE (used to compute the objective function)
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE FUSE_KGE_MODULE                           		   ! run model and compute the KGE
USE multiforce, only: ncid_forc                       ! NetCDF forcing file ID

IMPLICIT NONE
! input
INTEGER(I4B)                           :: NOPT        ! number of parameters
REAL(MSP), DIMENSION(100), INTENT(IN)  :: A            ! model parameter set - can be bumped up to 100 elements

! internal
REAL(SP), DIMENSION(:), ALLOCATABLE    :: SCE_PAR     ! sce parameter set
INTEGER(I4B)                           :: IERR        ! error code for allocate/deallocate
INTEGER(I4B)                           :: ERR         ! error code for fuse_kge
CHARACTER(LEN=256)                     :: MESSAGE     ! error message for fuse_kge
LOGICAL(LGT)                           :: OUTPUT_FLAG ! .TRUE. = write model time series
REAL(SP)                               :: KGE         ! KGE

! output
REAL(MSP)                              :: FUNCTN      ! objective function value

! ---------------------------------------------------------------------------------------
! get SCE parameter set
ALLOCATE(SCE_PAR(NOPT), STAT=IERR); IF (IERR.NE.0) STOP ' problem allocating space '
SCE_PAR(1:NOPT) = A(1:NOPT)  ! convert from MSP used in SCE to SP used in FUSE

OUTPUT_FLAG=.FALSE.   ! do not produce *runs.nc files only, param.nc files

CALL FUSE_KGE(SCE_PAR,.FALSE.,NCID_FORC,KGE,OUTPUT_FLAG,1) ! 2nd argument FALSE, always return KGE value

! deallocate parameter set
DEALLOCATE(SCE_PAR, STAT=IERR); IF (IERR.NE.0) STOP ' problem deallocating space '
print *, 'KGE =', KGE

! save objective function value
FUNCTN = -KGE
! ---------------------------------------------------------------------------------------
END FUNCTION FUNCTN
