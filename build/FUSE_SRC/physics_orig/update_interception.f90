SUBROUTINE UPDATE_INTERCEPTION(DT)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Cyril Thebault, 2026
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Updates the interception store over one forcing interval.
!
! PIN0      : precipitation entering the interception store [depth/time]
! EVAP_0    : evaporation from the interception store        [depth/time]
! PTHRU     : throughfall leaving the interception store      [depth/time]
! SINT_0    : interception storage                            [depth]
! MAXSINT_0 : maximum interception storage                    [depth]
! ---------------------------------------------------------------------------------------

USE nrtype
USE model_defn
USE model_defnames
USE multiparam
USE multiforce, ONLY: MFORCE
USE multistate
USE multi_flux

IMPLICIT NONE

REAL(WP), INTENT(IN) :: DT

REAL(WP) :: WATER_AVAILABLE
REAL(WP) :: WATER_REMAIN
REAL(WP) :: EVAP_AMOUNT
REAL(WP) :: PTHRU_AMOUNT

IF (DT <= 0._wp) THEN
  PRINT *, 'UPDATE_INTERCEPTION: DT must be greater than zero'
  STOP
END IF

SELECT CASE(SMODL%iINTRC)

 CASE(iopt_no_intrcep)

  M_FLUX%EVAP_0 = 0._wp
  M_FLUX%PTHRU  = M_FLUX%PIN0
  FSTATE%SINT_0 = 0._wp

 CASE(iopt_gr5h_intrc)

  IF (MPARAM%MAXSINT_0 < 0._wp) THEN
    PRINT *, 'UPDATE_INTERCEPTION: MAXSINT_0 must not be negative'
    STOP
  END IF

  WATER_AVAILABLE = MAX(0._wp, FSTATE%SINT_0)    &
                  + MAX(0._wp, M_FLUX%PIN0) * DT

  EVAP_AMOUNT = MIN(                              &
      MAX(0._wp, MFORCE%PET) * DT,               &
      WATER_AVAILABLE                            &
  )

  WATER_REMAIN = WATER_AVAILABLE - EVAP_AMOUNT

  PTHRU_AMOUNT = MAX(                             &
      0._wp,                                     &
      WATER_REMAIN - MPARAM%MAXSINT_0            &
  )

  FSTATE%SINT_0 = WATER_REMAIN - PTHRU_AMOUNT

  M_FLUX%EVAP_0 = EVAP_AMOUNT  / DT
  M_FLUX%PTHRU  = PTHRU_AMOUNT / DT

 CASE DEFAULT

  PRINT *, 'SMODL%iINTRC must be iopt_no_intrcep or iopt_gr5h_intrc'
  STOP

END SELECT

M_FLUX%EFF_PPT = M_FLUX%PTHRU

END SUBROUTINE UPDATE_INTERCEPTION