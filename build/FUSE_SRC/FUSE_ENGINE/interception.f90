 SUBROUTINE INTERCEPTION()
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Cyril Thebault to include interception, 7/2026
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Computes interception fluxes and updates the effective precipitation available to soil.
! ---------------------------------------------------------------------------------------
USE nrtype
USE model_defn
USE model_defnames
USE multiparam
USE multistate
USE multiforce
USE multi_flux

IMPLICIT NONE

! ---------------------------------------------------------------------------------------
SELECT CASE(SMODL%iINTRC)
 CASE(iopt_no_intrcep)
  M_FLUX%EVAP_0  = 0._SP
  M_FLUX%PTHRU = M_FLUX%PINC
 CASE(iopt_gr5h_intrc)
  M_FLUX%EVAP_0 = MIN(MFORCE%PET, M_FLUX%PINC + MSTATE%SINT_0/DELTIM)
  M_FLUX%PTHRU = MAX(0._SP, &
                     M_FLUX%PINC - (MPARAM%MAXSINT_0 - MSTATE%SINT_0)/DELTIM - M_FLUX%EVAP_0)
 CASE DEFAULT
  print *, "SMODL%iINTRC must be iopt_no_intrcep or iopt_gr5h_intrc"
  STOP
END SELECT

! keep existing FUSE convention: downstream modules consume EFF_PPT
M_FLUX%EFF_PPT = M_FLUX%PTHRU
! ---------------------------------------------------------------------------------------
END SUBROUTINE INTERCEPTION
