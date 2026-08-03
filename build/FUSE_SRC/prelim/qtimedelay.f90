SUBROUTINE QTIMEDELAY(info,err,message)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2007
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Computes the fraction of runoff in future time steps
! ---------------------------------------------------------------------------------------
! Modules Modified:
! -----------------
! MODULE multiparam -- runoff fractions stored in DPARAM%FRAC_FUTURE(:)
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE nr, ONLY : gammp                                  ! interface for the incomplete gamma function
USE model_defn                                        ! model definition structure
USE model_defnames
USE multiparam                                        ! model parameters
USE info_types, ONLY: fuse_info
USE multiroute, ONLY: FUTURE
IMPLICIT NONE
! dummies
integer(i4b),intent(out)::err
character(*),intent(out)::message
! locals
INTEGER(I4B)                           :: NTDH        ! maximum number of future time steps
REAL(WP)                               :: ALPHA       ! shape parameter
REAL(WP)                               :: ALAMB       ! scale parameter
INTEGER(I4B)                           :: JTIM        ! (loop through future time steps)
REAL(WP)                               :: TFUTURE     ! future time (units of days)
REAL(WP)                               :: CUMPROB     ! cumulative probability at JTIM
REAL(WP)                               :: PSAVE       ! cumulative probability at JTIM-1
TYPE(fuse_info), INTENT(IN)            :: info
INTEGER(I4B)                           :: ISTAT
! ---------------------------------------------------------------------------------------
err=0

! Compute the number of routing bins from the maximum routing horizon and the forcing timestep, both expressed in days.
NTDH = CEILING(TDH_MAX / info%time%deltim_days, KIND=I4B)

IF (NTDH < 2) THEN
  ERR = 100
  MESSAGE = 'f-QTIMEDELAY/at least two routing bins are required'
  RETURN
END IF

! Allocate or resize the runoff fractions.
IF (ALLOCATED(DPARAM%FRAC_FUTURE)) THEN
  IF (SIZE(DPARAM%FRAC_FUTURE) /= NTDH) THEN
    DEALLOCATE(DPARAM%FRAC_FUTURE, STAT=ISTAT)
    IF (ISTAT /= 0) THEN
      ERR = 100
      MESSAGE = 'f-QTIMEDELAY/cannot deallocate DPARAM%FRAC_FUTURE'
      RETURN
    END IF
  END IF
END IF

IF (.NOT.ALLOCATED(DPARAM%FRAC_FUTURE)) THEN
  ALLOCATE(DPARAM%FRAC_FUTURE(NTDH), STAT=ISTAT)
  IF (ISTAT /= 0) THEN
    ERR = 100
    MESSAGE = 'f-QTIMEDELAY/cannot allocate DPARAM%FRAC_FUTURE'
    RETURN
  END IF
END IF

! Allocate or resize the routing queue.
IF (ALLOCATED(FUTURE)) THEN
  IF (SIZE(FUTURE) /= NTDH) THEN
    DEALLOCATE(FUTURE, STAT=ISTAT)
    IF (ISTAT /= 0) THEN
      ERR = 100
      MESSAGE = 'f-QTIMEDELAY/cannot deallocate FUTURE'
      RETURN
    END IF
  END IF
END IF

IF (.NOT.ALLOCATED(FUTURE)) THEN
  ALLOCATE(FUTURE(NTDH), STAT=ISTAT)
  IF (ISTAT /= 0) THEN
    ERR = 100
    MESSAGE = 'f-QTIMEDELAY/cannot allocate FUTURE'
    RETURN
  END IF

  FUTURE = 0._WP
END IF

SELECT CASE(SMODL%iQ_TDH)
 CASE(iopt_rout_gamma) ! use a Gamma distribution with shape parameter = 2.5
  ALPHA = 2.5_WP                                             ! shape parameter

  !PRINT *, 'MPARAM= ', MPARAM

  ALAMB = ALPHA/MPARAM%TIMEDELAY                             ! scale parameter
  PSAVE = 0._WP                                              ! cumulative probability at JTIM-1
  NTDH = SIZE(DPARAM%FRAC_FUTURE)                            ! maximum number of future time steps
  ! loop through time steps and compute the fraction of runoff in future time steps
  DO JTIM=1,NTDH
   TFUTURE                   = REAL(JTIM, WP)*info%time%deltim_days          ! future time (units of days)
   CUMPROB                   = GAMMP(ALPHA,ALAMB*TFUTURE)    ! cumulative probability at JTIM
   DPARAM%FRAC_FUTURE(JTIM)  = MAX(0._WP, CUMPROB-PSAVE)     ! probability between JTIM-1 and JTIM
   PSAVE                     = CUMPROB                       ! cumulative probability at JTIM-1
   !WRITE(*,'(3(F11.5))') TFUTURE, DPARAM%FRAC_FUTURE(JTIM), CUMPROB
   IF(DPARAM%FRAC_FUTURE(JTIM)<EPSILON(1._WP))EXIT
  END DO
  DPARAM%NTDH_NEED = MIN(JTIM,NTDH)
  DPARAM%FRAC_FUTURE(DPARAM%NTDH_NEED+1:)=0._WP
  ! check there are enough bins
  IF (CUMPROB.LT.0.99_WP) THEN
   err=100; message='f-QTIMEDELAY/not enough bins in dparam%frac_future'
   return
  ENDIF
  ! ensure that the fractions sum to 1.0 (account for rounding errors, and not enough bins)
  DPARAM%FRAC_FUTURE(:) = DPARAM%FRAC_FUTURE(:) / SUM(DPARAM%FRAC_FUTURE(:))
 CASE(iopt_no_routing) ! no routing
  NTDH                       = SIZE(DPARAM%FRAC_FUTURE)
  DPARAM%NTDH_NEED           = 2
  DPARAM%FRAC_FUTURE(1)      = 1._WP
  DPARAM%FRAC_FUTURE(2:NTDH) = 0._WP
 CASE DEFAULT       ! check for errors
  err=100; message="f-QTIMEDELAY/SMODL%iQ_TDH must be either iopt_rout_gamma or iopt_no_routing"
  return
END SELECT
! ---------------------------------------------------------------------------------------
END SUBROUTINE QTIMEDELAY
