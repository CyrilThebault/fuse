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
  M_FLUX%PTHRU = M_FLUX%EFF_PPT
  M_FLUX%EINT  = 0._SP

 CASE(iopt_gr5h_intrc)
  ! temporary passive implementation
  M_FLUX%PTHRU = M_FLUX%EFF_PPT
  M_FLUX%EINT  = 0._SP

 CASE DEFAULT
  print *, "SMODL%iINTRC must be iopt_no_intrcep or iopt_gr5h_intrc"
  STOP
END SELECT

! keep existing FUSE convention: downstream modules consume EFF_PPT
M_FLUX%EFF_PPT = M_FLUX%PTHRU
! ---------------------------------------------------------------------------------------
END SUBROUTINE INTERCEPTION
