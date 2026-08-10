SUBROUTINE INITFLUXES(M_FLUX)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2007
! Modified by Brian Henn to include snow model, 6/2013
! Modified by Cyril Thébault to include interception, 7/2026
! Modified by Martyn Clark to transition to new data structures, 8/2026
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Set all fluxes to zero at the start of each time step
! ---------------------------------------------------------------------------------------
USE nrtype
use multi_flux_types,  only: FLUXES
implicit none
type(FLUXES), intent(out)  :: M_FLUX
! ---------------------------------------------------------------------------------------

M_FLUX%EVAP_0      = 0._wp
M_FLUX%PIN0        = 0._wp
M_FLUX%PTHRU       = 0._wp
M_FLUX%EFF_PPT     = 0._wp
M_FLUX%SATAREA     = 0._wp
M_FLUX%QSURF       = 0._wp
M_FLUX%EVAP_1A     = 0._wp
M_FLUX%EVAP_1B     = 0._wp
M_FLUX%EVAP_1      = 0._wp
M_FLUX%EVAP_2      = 0._wp
M_FLUX%RCHR2EXCS   = 0._wp
M_FLUX%TENS2FREE_1 = 0._wp
M_FLUX%TENS2FREE_2 = 0._wp
M_FLUX%QINTF_1     = 0._wp
M_FLUX%QPERC_12    = 0._wp
M_FLUX%QBASE_2     = 0._wp
M_FLUX%QBASE_2A    = 0._wp
M_FLUX%QBASE_2B    = 0._wp
M_FLUX%OFLOW_1     = 0._wp
M_FLUX%OFLOW_2     = 0._wp
M_FLUX%OFLOW_2A    = 0._wp
M_FLUX%OFLOW_2B    = 0._wp

M_FLUX%ERR_WATR_1  = 0._wp
M_FLUX%ERR_TENS_1  = 0._wp
M_FLUX%ERR_FREE_1  = 0._wp
M_FLUX%ERR_TENS_1A = 0._wp
M_FLUX%ERR_TENS_1B = 0._wp
M_FLUX%ERR_WATR_2  = 0._wp
M_FLUX%ERR_TENS_2  = 0._wp
M_FLUX%ERR_FREE_2  = 0._wp
M_FLUX%ERR_FREE_2A = 0._wp
M_FLUX%ERR_FREE_2B = 0._wp
M_FLUX%CHK_TIME    = 0._wp
! ---------------------------------------------------------------------------------------
END SUBROUTINE INITFLUXES
