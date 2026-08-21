SUBROUTINE INIT_STATS(work)
! ----------------------------------------------------------------------------------------
! Creator:
!   Martyn Clark, 2009
!
! ----------------------------------------------------------------------------------------
! Purpose:
!   Used to initialize summary statistics
!
! ----------------------------------------------------------------------------------------
! Future revisions:
!
!   (add other summary statistics)
!
! ----------------------------------------------------------------------------------------
USE work_types, only: fuse_work
USE model_numerix, only: PRB_NSUBS
IMPLICIT NONE
type(fuse_work), intent(inout)    :: work
! ----------------------------------------------------------------------------------------
! initialize numerical statistics
work%run%stats%NUM_FUNCS     = 0 
work%run%stats%NUM_JACOBIAN  = 0
work%run%stats%NUMSUB_ACCEPT = 0
work%run%stats%NUMSUB_REJECT = 0
work%run%stats%NUMSUB_NOCONV = 0
! initialize probability distributions
PRB_NSUBS(:) = 0
! ----------------------------------------------------------------------------------------
END SUBROUTINE INIT_STATS
