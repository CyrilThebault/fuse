FUNCTION FUNCTN(NOPT,A)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2009
! Modified by Cyril Thébault to allow different metrics as objective function, 2024
! Modified by Cyril Thébault to allow different transformations on parameters, 2026
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Wrapper for SCE (used to compute the objective function)
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE fuse_metric_module                                ! run model and compute the metric chosen as objective function
USE multiforce, only: ncid_forc                       ! NetCDF forcing file ID
USE fuse_fileManager,only:METRIC, TRANSFO             ! metric and transformation requested in the filemanager

! parameter transformation utilities and calibration state
USE parameter_transform_module, ONLY: &               ! convert parameters back to physical space
  vector_to_physical_space
USE calibration_data_module, ONLY: &                  ! transformation codes for current calibration
  CALIB_TRANSFORM_CODES

IMPLICIT NONE
! input
INTEGER(I4B)                           :: NOPT        ! number of parameters
REAL(MSP), DIMENSION(100), INTENT(IN)  :: A            ! model parameter set - can be bumped up to 100 elements

! internal
REAL(SP), DIMENSION(:), ALLOCATABLE    :: SCE_PAR     ! sce parameter set
INTEGER(I4B)                           :: IERR        ! error code for allocate/deallocate
INTEGER(I4B)                           :: ERR         ! error code for fuse_metric
CHARACTER(LEN=256)                     :: MESSAGE     ! error message for fuse_metric
LOGICAL(LGT)                           :: OUTPUT_FLAG ! .TRUE. = write model time series
REAL(SP)                               :: METRIC_VAL  ! value of the metric chosen as objective function

! parameter transformation
REAL(MSP), DIMENSION(:), ALLOCATABLE   :: PHYSICAL_PAR_MSP ! parameters in physical space

! output
REAL(MSP)                              :: FUNCTN      ! objective function value

! ---------------------------------------------------------------------------------------
! allocate parameter arrays
ALLOCATE(SCE_PAR(NOPT), PHYSICAL_PAR_MSP(NOPT), STAT=IERR)
IF (IERR.NE.0) STOP ' problem allocating parameter arrays '

! ensure that transformation codes were initialized by the calibration driver
IF (.NOT.ALLOCATED(CALIB_TRANSFORM_CODES)) THEN
  STOP ' calibration transformation codes have not been initialized '
END IF

! ensure consistency between SCE parameters and calibration metadata
IF (SIZE(CALIB_TRANSFORM_CODES).NE.NOPT) THEN
  STOP ' inconsistent number of calibration transformation codes '
END IF

! convert parameters from the SCE search space to the physical FUSE space
CALL vector_to_physical_space( &
     A(1:NOPT),                &
     CALIB_TRANSFORM_CODES,    &
     PHYSICAL_PAR_MSP,         &
     ERR,                      &
     MESSAGE)

IF (ERR.NE.0) THEN
  WRITE(*,'(A)') TRIM(MESSAGE)
  STOP ' unable to transform parameters to physical space '
END IF

! convert from MSP used in SCE to SP used in FUSE
SCE_PAR(1:NOPT) = REAL(PHYSICAL_PAR_MSP(1:NOPT), SP)

OUTPUT_FLAG=.FALSE.   ! do not produce *runs.nc files only, param.nc files

CALL FUSE_METRIC(SCE_PAR,.FALSE.,NCID_FORC,METRIC_VAL,OUTPUT_FLAG,1) ! 2nd argument FALSE, always return METRIC value

! deallocate parameter arrays
DEALLOCATE(SCE_PAR, PHYSICAL_PAR_MSP, STAT=IERR)
IF (IERR.NE.0) STOP ' problem deallocating parameter arrays '

! save objective function value: SCE is a minimization algorithm
IF (METRIC=="KGE" .OR. METRIC=="KGEP" .OR. METRIC=="NSE") THEN
  FUNCTN = -METRIC_VAL
ELSE IF (METRIC=="MAE" .OR. METRIC=="RMSE" ) THEN
  FUNCTN = METRIC_VAL
ELSE 
   STOP 'The requested metric is not available in metrics module'
END IF

! ---------------------------------------------------------------------------------------
END FUNCTION FUNCTN
