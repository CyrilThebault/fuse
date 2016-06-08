SUBROUTINE INIT_STATE(FRAC)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2007
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Initialize model states at fraction (FRAC) of capacity
! ---------------------------------------------------------------------------------------
! Modules Modified:
! -----------------
! Model states in MODULE multistate
! ---------------------------------------------------------------------------------------
USE multiparam                                        ! model parameters
USE multistate                                        ! model states
USE multiroute                                        ! routed runoff
IMPLICIT NONE
REAL(SP), INTENT(IN)                  :: FRAC         ! fraction of capacity
! ---------------------------------------------------------------------------------------
! (upper layer)
FSTATE%TENS_1A = DPARAM%MAXTENS_1A * FRAC
FSTATE%TENS_1B = DPARAM%MAXTENS_1B * FRAC
FSTATE%TENS_1  = DPARAM%MAXTENS_1  * FRAC
FSTATE%FREE_1  = DPARAM%MAXFREE_1  * FRAC
FSTATE%WATR_1  = MPARAM%MAXWATR_1  * FRAC
! (lower layer)
FSTATE%TENS_2  = DPARAM%MAXTENS_2  * FRAC
FSTATE%FREE_2  = DPARAM%MAXFREE_2  * FRAC
FSTATE%FREE_2A = DPARAM%MAXFREE_2A * FRAC
FSTATE%FREE_2B = DPARAM%MAXFREE_2B * FRAC
FSTATE%WATR_2  = MPARAM%MAXWATR_2  * FRAC
! (routed runoff)
FUTURE         = 0._sp
! ---------------------------------------------------------------------------------------
END SUBROUTINE INIT_STATE
