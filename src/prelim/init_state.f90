module init_state_module

  use nrtype
  
  use multiparam_types, only: PARADJ
  use multiparam_types, only: PARDVD
  use multistate_types, only: STATEV
  use multibands_types, only: BANDS_VAR

  implicit none

  private
  public :: INIT_STATE

contains

  SUBROUTINE INIT_STATE(FRAC, MPARAM, DPARAM, FSTATE)
  ! ---------------------------------------------------------------------------------------
  ! Creator:
  ! --------
  ! Martyn Clark, 2007
  ! Modified by Brian Henn to include snow model, 6/2013
  ! Modified by Cyril Thebault to include interception, 7/2026
  ! Modified by Martyn Clark to use new data structures, 8/2026
  ! ---------------------------------------------------------------------------------------
  ! Purpose:
  ! --------
  ! Initialize model states at fraction (FRAC) of capacity
  ! ---------------------------------------------------------------------------------------
  ! Modules Modified:
  ! -----------------
  ! Model states in MODULE multistate
  ! ---------------------------------------------------------------------------------------
  IMPLICIT NONE
  REAL(WP)         , INTENT(IN)      :: FRAC        ! fraction of capacity
  type(PARADJ)     , intent(in)      :: MPARAM      ! adjustable model parameters (time delay)
  type(PARDVD)     , intent(in)      :: DPARAM      ! derived model parameters (FRAC_FUTURE, )
  type(STATEV)     , intent(inout)   :: FSTATE      ! model states
  ! ---------------------------------------------------------------------------------------
  ! interception layer
  FSTATE%SINT_0 = 0._wp
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
  ! ---------------------------------------------------------------------------------------
  END SUBROUTINE INIT_STATE

end module init_state_module
