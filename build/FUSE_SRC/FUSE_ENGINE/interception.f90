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
REAL(SP) :: DT
REAL(SP) :: P_STEP
REAL(SP) :: PET_STEP
REAL(SP) :: EINT_STEP
REAL(SP) :: PTHRU_STEP
! ---------------------------------------------------------------------------------------
SELECT CASE(SMODL%iINTRC)
 CASE(iopt_no_intrcep)
  M_FLUX%EVAP_0  = 0._SP
  M_FLUX%PTHRU = M_FLUX%PINC
 CASE(iopt_gr5h_intrc)
  DT = HSTATE%STEP

  P_STEP   = M_FLUX%PINC * DT
  PET_STEP = MFORCE%PET  * DT

  EINT_STEP = MIN(PET_STEP, P_STEP + MSTATE%SINT_0)

  PTHRU_STEP = MAX(0._SP, &
                   P_STEP - (MPARAM%MAXSINT_0 - MSTATE%SINT_0) - EINT_STEP)                

  M_FLUX%EVAP_0 = EINT_STEP / DT
  M_FLUX%PTHRU = PTHRU_STEP / DT

 CASE DEFAULT
  print *, "SMODL%iINTRC must be iopt_no_intrcep or iopt_gr5h_intrc"
  STOP
END SELECT

! keep existing FUSE convention: downstream modules consume EFF_PPT
M_FLUX%EFF_PPT = M_FLUX%PTHRU
! ---------------------------------------------------------------------------------------
END SUBROUTINE INTERCEPTION
