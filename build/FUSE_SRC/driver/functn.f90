FUNCTION FUNCTN(NOPT,A)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2009
! Modified by Cyril Thébault to allow different metrics as objective function, 2024
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Wrapper for SCE (used to compute the objective function)
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE sce_callback_context, only: ctx                   ! access FUSE data structures
USE fuse_evaluate_module, only: fuse_evaluate         ! run model and compute the metric chosen as objective function
USE multiforce, only: ncid_forc                       ! NetCDF forcing file ID
USE fuse_fileManager,only:METRIC, TRANSFO             ! metric and transformation requested in the filemanager
USE fuse_globaldata, only: nFUSE_eval                      ! # fuse evaluations

IMPLICIT NONE
! input
INTEGER(I4B)                           :: NOPT        ! number of parameters
REAL(MSP), DIMENSION(100), INTENT(IN)  :: A            ! model parameter set - can be bumped up to 100 elements

! internal
REAL(SP), DIMENSION(NOPT)              :: SCE_PAR     ! sce parameter set
INTEGER(I4B)                           :: IERR        ! error code for allocate/deallocate
INTEGER(I4B)                           :: ERR         ! error code for fuse_metric
CHARACTER(LEN=256)                     :: MESSAGE     ! error message for fuse_metric
LOGICAL(LGT)                           :: OUTPUT_FLAG ! .TRUE. = write model time series
REAL(SP)                               :: METRIC_VAL  ! value of the metric chosen as objective function

! output
REAL(MSP)                              :: FUNCTN      ! objective function value

! ---------------------------------------------------------------------------------------

nFUSE_eval = nFUSE_eval + 1

! get SCE parameter set
SCE_PAR(1:NOPT) = A(1:NOPT)  ! convert from MSP used in SCE to SP used in FUSE
OUTPUT_FLAG=.FALSE.          ! do not produce *runs.nc files only, param.nc files

CALL FUSE_evaluate(SCE_PAR, ctx%info, ctx%work, ctx%domain, OUTPUT_FLAG, METRIC_VAL)

! save objective function value: SCE is a minimization algorithm
select case(metric)
 case ("KGE", "KGEP", "NSE"); FUNCTN = -METRIC_VAL
 case ("MAE", "RMSE");        FUNCTN = METRIC_VAL
 case default
  STOP 'The requested metric is not available in metrics module'
end select

! ---------------------------------------------------------------------------------------
END FUNCTION FUNCTN
