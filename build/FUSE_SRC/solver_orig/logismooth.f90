PURE FUNCTION LOGISMOOTH(STATE,STATE_MAX,PSMOOTH)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2007
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Uses a logistic function to smooth the threshold at the top of a bucket
! ---------------------------------------------------------------------------------------
USE nrtype
IMPLICIT NONE
REAL(WP), INTENT(IN)                   :: STATE       ! model state
REAL(WP), INTENT(IN)                   :: STATE_MAX   ! maximum model state
REAL(WP), INTENT(IN)                   :: PSMOOTH     ! smoothing parameter (fraction of state)
REAL(WP)                               :: ASMOOTH     ! actual smoothing
REAL(WP)                               :: LOGISMOOTH  ! FUNCTION name
! ---------------------------------------------------------------------------------------
ASMOOTH = PSMOOTH*STATE_MAX                           ! actual smoothing
LOGISMOOTH = 1._WP / ( 1._WP + EXP(-(STATE - (STATE_MAX - ASMOOTH*5._WP) ) / ASMOOTH) )
! ---------------------------------------------------------------------------------------
END FUNCTION LOGISMOOTH
