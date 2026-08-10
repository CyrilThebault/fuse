module QTIMEDELAY_module

  implicit none

  private
  public :: QTIMEDELAY

contains

  SUBROUTINE QTIMEDELAY(info, MPARAM, DPARAM, FUTURE, ierr, message)
  ! ---------------------------------------------------------------------------------------
  ! Creator:
  ! --------
  ! Martyn Clark, 2007, modified 08/2026 to use new data structures
  ! ---------------------------------------------------------------------------------------
  ! Purpose:
  ! --------
  ! Computes the fraction of runoff in future time steps
  ! ---------------------------------------------------------------------------------------
  
  ! data types
  use nrtype                                            ! variable types, etc.
  use info_types,       only: fuse_info
  use multiparam_types, only: PARADJ            
  use multiparam_types, only: PARDVD
  
  ! model options
  USE model_defn                                        ! model definition structure
  USE model_defnames
  
  ! recipes
  USE nr, ONLY : gammp                                  ! interface for the incomplete gamma function
  
  IMPLICIT NONE
  
  ! dummies
  type(fuse_info)      , intent(in)     :: info
  type(PARADJ)         , intent(in)     :: MPARAM      ! adjustable model parameters (time delay) 
  type(PARDVD)         , intent(inout)  :: DPARAM      ! derived model parameters (FRAC_FUTURE, )
  real(wp), allocatable, intent(inout)  :: FUTURE(:)   ! rolling routing convolution queue
  
  integer(i4b)     , intent(out)        :: ierr
  character(*)     , intent(out)        :: message
  
  ! locals
  INTEGER(I4B)                          :: NTDH        ! maximum number of future time steps
  REAL(WP)                              :: ALPHA       ! shape parameter
  REAL(WP)                              :: ALAMB       ! scale parameter
  INTEGER(I4B)                          :: JTIM        ! (loop through future time steps)
  REAL(WP)                              :: TFUTURE     ! future time (units of days)
  REAL(WP)                              :: CUMPROB     ! cumulative probability at JTIM
  REAL(WP)                              :: PSAVE       ! cumulative probability at JTIM-1
  INTEGER(I4B)                          :: ISTAT
  ! ---------------------------------------------------------------------------------------
  ierr = 0
  message = "QTIMEDELAY/"
  
  ! Compute the number of routing bins from the maximum routing horizon and the forcing timestep, both expressed in days.
  NTDH = CEILING(TDH_MAX / info%time%deltim_days, KIND=I4B)
  
  IF (NTDH < 2) THEN
    MESSAGE = trim(message)//'at least two routing bins are required'
    IERR = 100; RETURN
  END IF
  
  ! Allocate or resize runoff fractions
  if (allocated(DPARAM%FRAC_FUTURE)) then
    if (size(DPARAM%FRAC_FUTURE) /= NTDH) deallocate(DPARAM%FRAC_FUTURE)
  end if
  
  if (.not. allocated(DPARAM%FRAC_FUTURE)) then
    allocate(DPARAM%FRAC_FUTURE(NTDH), stat=istat)
    if (istat /= 0) then
      message = trim(message)//'cannot allocate DPARAM%FRAC_FUTURE'
      ierr = 100; return
    end if
  end if
  
  ! Allocate or resize routing queue
  if (allocated(FUTURE)) then
    if (size(FUTURE) /= NTDH) deallocate(FUTURE)
  end if
  
  if (.not. allocated(FUTURE)) then
  
    allocate(FUTURE(NTDH), stat=istat)
    if (istat /= 0) then
      message = trim(message)//'cannot allocate FUTURE'
      ierr = 100; return
    end if
  
    FUTURE = 0._wp
  
  end if
  
  SELECT CASE(SMODL%iQ_TDH)
  
    CASE (iopt_rout_gamma) ! use a Gamma distribution with shape parameter = 2.5
    
      ALPHA = 2.5_WP                                             ! shape parameter
      ALAMB = ALPHA/MPARAM%TIMEDELAY                             ! scale parameter
      NTDH  = SIZE(DPARAM%FRAC_FUTURE)                           ! maximum number of future time steps
      
      PSAVE = 0._wp
     
      do JTIM = 1, NTDH
     
        TFUTURE = real(JTIM,wp) * info%time%deltim_days
        CUMPROB = GAMMP(ALPHA, ALAMB*TFUTURE)
     
        DPARAM%FRAC_FUTURE(JTIM) = max(0._wp, CUMPROB - PSAVE)
        PSAVE = CUMPROB
        
        !WRITE(*,'(3(F11.5))') TFUTURE, DPARAM%FRAC_FUTURE(JTIM), CUMPROB 
        
        if (CUMPROB >= 0.999_wp) exit
     
      end do
     
      DPARAM%NTDH_NEED = JTIM
     
      ! make sure enough probability mass was captured
      if (CUMPROB < 0.99_wp) then
        message = trim(message)//'not enough bins in dparam%frac_future'
        ierr = 100; return
      end if
     
      ! renormalize retained routing fractions
      DPARAM%FRAC_FUTURE(1:DPARAM%NTDH_NEED) = DPARAM%FRAC_FUTURE(1:DPARAM%NTDH_NEED) / &
                                           sum(DPARAM%FRAC_FUTURE(1:DPARAM%NTDH_NEED))
     
    CASE (iopt_no_routing) ! no routing
    
      NTDH                       = SIZE(DPARAM%FRAC_FUTURE)
      DPARAM%NTDH_NEED           = 2
      DPARAM%FRAC_FUTURE(1)      = 1._WP
      DPARAM%FRAC_FUTURE(2:NTDH) = 0._WP
   
    CASE DEFAULT       ! check for nerrors
      message=trim(message)//"SMODL%iQ_TDH must be either iopt_rout_gamma or iopt_no_routing"
      ierr=100; return
  
  END SELECT
  ! ---------------------------------------------------------------------------------------
  END SUBROUTINE QTIMEDELAY

end module QTIMEDELAY_module
