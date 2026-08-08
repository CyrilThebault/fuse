SUBROUTINE COMP_STATS(work)
! ----------------------------------------------------------------------------------------
! Creator:
!   Martyn Clark, 2009
!
! ----------------------------------------------------------------------------------------
! Purpose:
!   Used to compute summary statistics of model output
!
! ----------------------------------------------------------------------------------------
! Future revisions:
!
!   (add other summary statistics)
!
! ----------------------------------------------------------------------------------------
Use nrtype
USE work_types, only: fuse_work
USE model_numerix, only: NUM_FUNCS, NUM_JACOBIAN, NUMSUB_ACCEPT, NUMSUB_REJECT, NUMSUB_NOCONV
USE model_numerix, only: MAXNUM_ITERNS, ORD_NSUBS, PRB_NSUBS
IMPLICIT NONE
type(fuse_work), intent(inout)    :: work  ! structures that depend on nState/nPar
! ----------------------------------------------------------------------------------------
! compute numerical stats
work%run%stats%NUM_FUNCS     = work%run%stats%NUM_FUNCS     + REAL(NUM_FUNCS, KIND(WP))     ! number of function calls
work%run%stats%NUM_JACOBIAN  = work%run%stats%NUM_JACOBIAN  + REAL(NUM_JACOBIAN, KIND(WP))  ! number of times Jacobian is calculated
work%run%stats%NUMSUB_ACCEPT = work%run%stats%NUMSUB_ACCEPT + REAL(NUMSUB_ACCEPT, KIND(WP)) ! number of sub-steps accepted (taken)
work%run%stats%NUMSUB_REJECT = work%run%stats%NUMSUB_REJECT + REAL(NUMSUB_REJECT, KIND(WP)) ! number of sub-steps tried but rejected
work%run%stats%NUMSUB_NOCONV = work%run%stats%NUMSUB_NOCONV + REAL(NUMSUB_NOCONV, KIND(WP)) ! number of sub-steps tried that did not converge
! compute maximum number of iterations
IF (MAXNUM_ITERNS > work%run%stats%MAXNUM_ITERNS) work%run%stats%MAXNUM_ITERNS = MAXNUM_ITERNS
! compute probability distributions
WHERE(ORD_NSUBS.GE.NUMSUB_ACCEPT) PRB_NSUBS = PRB_NSUBS + 1
! ----------------------------------------------------------------------------------------
END SUBROUTINE COMP_STATS
