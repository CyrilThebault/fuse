! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Cyril Thébault, 2026
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Store information associated with the current calibration.
!
! This module separates the calibration state from the mathematical parameter
! transformation utilities defined in parameter_transfo.f90.
! ---------------------------------------------------------------------------------------

MODULE calibration_data_module

  USE nrtype, ONLY: I4B

  IMPLICIT NONE
  PRIVATE

  ! transformation code associated with each calibrated parameter
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE, PUBLIC :: CALIB_TRANSFORM_CODES

  PUBLIC :: INIT_CALIBRATION_TRANSFORMS
  PUBLIC :: CLEAR_CALIBRATION_TRANSFORMS

CONTAINS


  ! -------------------------------------------------------------------------------------
  ! Initialize the transformation codes used during the current calibration.
  ! -------------------------------------------------------------------------------------
  SUBROUTINE INIT_CALIBRATION_TRANSFORMS(NPAR)

    INTEGER(I4B), INTENT(IN) :: NPAR
    INTEGER(I4B)             :: IERR

    ! remove transformation codes from a previous calibration, if necessary
    IF (ALLOCATED(CALIB_TRANSFORM_CODES)) THEN

      DEALLOCATE(CALIB_TRANSFORM_CODES, STAT=IERR)

      IF (IERR.NE.0) THEN
        STOP ' unable to deallocate calibration transformation codes '
      END IF

    END IF

    ! check that the requested number of parameters is valid
    IF (NPAR.LE.0) THEN
      STOP ' invalid number of calibration transformation codes '
    END IF

    ! allocate transformation codes for the current calibration
    ALLOCATE(CALIB_TRANSFORM_CODES(NPAR), STAT=IERR)

    IF (IERR.NE.0) THEN
      STOP ' unable to allocate calibration transformation codes '
    END IF

    ! initialize all parameters with the identity transformation
    CALIB_TRANSFORM_CODES(:) = 0_I4B

  END SUBROUTINE INIT_CALIBRATION_TRANSFORMS
  

  ! -------------------------------------------------------------------------------------
  ! Clear the transformation codes after calibration.
  ! -------------------------------------------------------------------------------------
  SUBROUTINE CLEAR_CALIBRATION_TRANSFORMS()

    INTEGER(I4B) :: IERR

    IF (ALLOCATED(CALIB_TRANSFORM_CODES)) THEN

      DEALLOCATE(CALIB_TRANSFORM_CODES, STAT=IERR)

      IF (IERR.NE.0) THEN
        STOP ' unable to clear calibration transformation codes '
      END IF

    END IF

  END SUBROUTINE CLEAR_CALIBRATION_TRANSFORMS

END MODULE calibration_data_module
