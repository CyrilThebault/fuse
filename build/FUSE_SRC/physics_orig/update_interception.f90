SUBROUTINE UPDATE_INTERCEPTION(DT, IERR, MESSAGE)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Cyril Thebault, 2026
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Updates the interception store.
!
! PIN0      : precipitation entering the interception store   [depth/time]
! EVAP_0    : evaporation from the interception store         [depth/time]
! PTHRU     : throughfall leaving the interception store      [depth/time]
! SINT_0    : interception storage                            [depth]
! REFSINT_0 : characteristic interception storage             [depth]
! ---------------------------------------------------------------------------------------

USE nrtype
USE model_defn
USE model_defnames
USE multiparam
USE multiforce, ONLY: MFORCE
USE multistate
USE multi_flux
USE model_numerix, ONLY: ERR_ITER_FUNC, ERR_ITER_DX, NITER_TOTAL

IMPLICIT NONE

! Input
REAL(WP), INTENT(IN) :: DT

! Output
INTEGER(I4B), INTENT(OUT) :: IERR
CHARACTER(*), INTENT(OUT) :: MESSAGE

! Local variables
INTEGER(I4B) :: ITER

REAL(WP) :: S_OLD
REAL(WP) :: S_NEW
REAL(WP) :: S_TRIAL
REAL(WP) :: S_LO
REAL(WP) :: S_HI

REAL(WP) :: PRECIP
REAL(WP) :: PET

REAL(WP), PARAMETER :: K_SMOOTH = 1.e-5_wp    ! It may need to be flexible based on REFSINT_0; e.g. K_SMOOTH=MAX(1.e-6_wp, 0.01_wp * MPARAM%REFSINT_0)
REAL(WP), PARAMETER :: EPS_WET  = 1.e-6_wp

REAL(WP) :: X_LOGISTIC
REAL(WP) :: EXP_X
REAL(WP) :: PHI
REAL(WP) :: WET_FRAC

REAL(WP) :: DPHI_DS
REAL(WP) :: DWET_DS

REAL(WP) :: RESIDUAL
REAL(WP) :: DRESIDUAL
REAL(WP) :: NEWTON_STEP

LOGICAL(LGT) :: CONVERGED

! ---------------------------------------------------------------------------------------
! Initialise error handling
! ---------------------------------------------------------------------------------------
IERR = 0
MESSAGE = 'UPDATE_INTERCEPTION/'
CONVERGED = .FALSE.

IF (DT <= 0._wp) THEN
  IERR = 1
  MESSAGE = TRIM(MESSAGE)//'DT must be greater than zero'
  RETURN
END IF

! ---------------------------------------------------------------------------------------
! Interception options
! ---------------------------------------------------------------------------------------
SELECT CASE(SMODL%iINTRC)

 CASE(iopt_no_intrcep)

  FSTATE%SINT_0 = 0._wp

  M_FLUX%EVAP_0 = 0._wp
  M_FLUX%PTHRU  = M_FLUX%PIN0
  M_FLUX%EFF_PPT = M_FLUX%PTHRU

  RETURN

 CASE(iopt_gr5h_intrc)

  IF (MPARAM%REFSINT_0 < 0._wp) THEN
    IERR = 1
    MESSAGE = TRIM(MESSAGE)//'REFSINT_0 must not be negative'
    RETURN
  END IF

 CASE DEFAULT

  IERR = 1
  MESSAGE = TRIM(MESSAGE)//'SMODL%iINTRC must be iopt_no_intrcep or iopt_gr5h_intrc'
  RETURN

END SELECT

! ---------------------------------------------------------------------------------------
! Define forcing and smoothing scales
! ---------------------------------------------------------------------------------------

PRECIP = MAX(0._wp, M_FLUX%PIN0)
PET    = MAX(0._wp, MFORCE%PET)

S_OLD = MAX(0._wp, FSTATE%SINT_0)

! The storage cannot become negative and cannot exceed the initial storage
! plus all incoming precipitation when evaporation and throughfall are omitted.
S_LO = 0._wp
S_HI = S_OLD + PRECIP * DT

! Handle the trivial dry case explicitly.
IF (S_HI <= ERR_ITER_DX) THEN
  S_NEW = 0._wp
  CONVERGED = .TRUE.
ELSE
  ! Start Newton from the old storage, constrained to the bracket.
  S_NEW = MIN(MAX(S_OLD, S_LO), S_HI)
END IF

! ---------------------------------------------------------------------------------------
! Solve the backward-Euler residual
!
! R(S) = S - S_OLD - DT * [PRECIP*(1-PHI(S)) - PET*WET_FRAC(S)]
! ---------------------------------------------------------------------------------------
IF (.NOT.CONVERGED) THEN

  DO ITER = 1, NITER_TOTAL

    ! Stable evaluation of the logistic function:
    ! PHI = 1 / (1 + exp(-(S-REFSINT_0)/K_SMOOTH))
    X_LOGISTIC = (S_NEW - MPARAM%REFSINT_0) / K_SMOOTH

    IF (X_LOGISTIC >= 0._wp) THEN
      PHI = 1._wp / (1._wp + EXP(-X_LOGISTIC))
    ELSE
      EXP_X = EXP(X_LOGISTIC)
      PHI = EXP_X / (1._wp + EXP_X)
    END IF

    WET_FRAC = S_NEW / (S_NEW + EPS_WET)

    DPHI_DS = PHI * (1._wp - PHI) / K_SMOOTH
    DWET_DS = EPS_WET / (S_NEW + EPS_WET)**2

    RESIDUAL = S_NEW - S_OLD - DT * (PRECIP * (1._wp - PHI) - PET * WET_FRAC)

    DRESIDUAL = 1._wp + DT * (PRECIP * DPHI_DS + PET * DWET_DS)

    IF (ABS(RESIDUAL) <= ERR_ITER_FUNC) THEN
      CONVERGED = .TRUE.
      EXIT
    END IF

    ! The residual is monotonic. Update the root bracket first.
    IF (RESIDUAL < 0._wp) THEN
      S_LO = S_NEW
    ELSE
      S_HI = S_NEW
    END IF

    ! Newton candidate.
    NEWTON_STEP = -RESIDUAL / DRESIDUAL
    S_TRIAL = S_NEW + NEWTON_STEP

    ! Safeguard Newton with bisection.
    IF (S_TRIAL <= S_LO .OR. S_TRIAL >= S_HI) THEN
      S_TRIAL = 0.5_wp * (S_LO + S_HI)
    END IF

    IF (ABS(S_TRIAL - S_NEW) <= ERR_ITER_DX) THEN
      S_NEW = S_TRIAL
      CONVERGED = .TRUE.
      EXIT
    END IF

    S_NEW = S_TRIAL

  END DO

END IF

IF (.NOT.CONVERGED) THEN
  IERR = 1
  MESSAGE = TRIM(MESSAGE)// &
            'Newton-bisection failed to converge'
  RETURN
END IF

! ---------------------------------------------------------------------------------------
! Evaluate fluxes at the converged end-of-step storage
! ---------------------------------------------------------------------------------------
X_LOGISTIC = (S_NEW - MPARAM%REFSINT_0) / K_SMOOTH

IF (X_LOGISTIC >= 0._wp) THEN
  PHI = 1._wp / (1._wp + EXP(-X_LOGISTIC))
ELSE
  EXP_X = EXP(X_LOGISTIC)
  PHI = EXP_X / (1._wp + EXP_X)
END IF

WET_FRAC = S_NEW / (S_NEW + EPS_WET)

FSTATE%SINT_0 = MAX(0._wp, S_NEW)

M_FLUX%PTHRU  = PRECIP * PHI
M_FLUX%EVAP_0 = PET * WET_FRAC
M_FLUX%EFF_PPT = M_FLUX%PTHRU


END SUBROUTINE UPDATE_INTERCEPTION