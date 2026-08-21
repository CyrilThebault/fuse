module Q_OVERLAND_module

  use nrtype
  

  implicit none

  private
  public :: Q_OVERLAND

contains

  SUBROUTINE Q_OVERLAND(DPARAM, &
                        W_FLUX, &
                        FUTURE, &
                        MROUTE, &
                        ierr, message)
  ! ---------------------------------------------------------------------------------------
  ! Creator:
  ! --------
  ! Martyn Clark, 2007
  ! --------
  ! History
  ! 5 June 2013 AD: Modified by David McInerney to merge array loop operations
  ! 5 June 2013 AD: Modified by Dmitri Kavetski to avoid zero-element operations
  ! ---------------------------------------------------------------------------------------
  ! Purpose:
  ! --------
  ! Computes the time delay in runoff in a basin (places runoff in future time steps)
  ! ---------------------------------------------------------------------------------------
  ! Modules Modified:
  ! -----------------
  ! MODULE multiroute -- places runoff in array FUTURE(:)RUNOFF
  ! ---------------------------------------------------------------------------------------
  
  ! data types
  USE multiparam_types,  only: PARDVD
  USE multi_flux_types,  only: FLUXES
  USE multiroute_types,  only: RUNOFF
  
  ! model options
  use model_defn,   only: SMODL
  use model_defnames
  
  IMPLICIT NONE
  
  type(PARDVD)         , intent(in)     :: DPARAM      ! derived model parameters (FRAC_FUTURE)
  type(FLUXES)         , intent(in)     :: W_FLUX      ! flux dats structure (contributions to runoff)
  real(wp)             , intent(inout)  :: FUTURE(:)   ! rolling routing convolution queue
  type(RUNOFF)         , intent(out)    :: MROUTE      ! routing data structure

  integer(i4b)     , intent(out)        :: ierr
  character(len=*) , intent(out)        :: message
  
  ! locals
  INTEGER(I4B)                          :: NTDH        ! maximum number of future time steps
  INTEGER(I4B)                          :: MTDH        ! reduced maximum number of future time steps
  INTEGER(I4B)                          :: JTIM        ! (loop through future time steps)
  REAL(WP), PARAMETER                   :: SNEG=-1.e-5 ! small negative number, used for checking
  LOGICAL, PARAMETER                    :: USE_NTDH_NEED=.TRUE. ! flag to use NTDH_NEED to reduce array operations (loop length)
  
  ierr = 0
  message = 'Q_OVERLAND/'
  
  ! ---------------------------------------------------------------------------------------
  ! compute total runoff (sum of surface runoff, overflow, interflow, and baseflow)
  MROUTE%Q_INSTNT = ( W_FLUX%QSURF   + &
                      W_FLUX%OFLOW_1 + &
                      W_FLUX%QINTF_1 + &
                      W_FLUX%OFLOW_2 + &
                      W_FLUX%QBASE_2 )
  ! ---------------------------------------------------------------------------------------
  SELECT CASE(SMODL%iQ_TDH)
  
   CASE (iopt_rout_gamma)             ! use a Gamma distribution with shape parameter = 2.5
    
    NTDH = SIZE(DPARAM%FRAC_FUTURE)  ! maximum number of future time steps
    MTDH = MERGE(DPARAM%NTDH_NEED, NTDH, USE_NTDH_NEED) ! reduce array operations (loop length)
    
    ! routed flow
    MROUTE%Q_ROUTED = FUTURE(1) + MROUTE%Q_INSTNT * DPARAM%FRAC_FUTURE(1)
    
    ! place runoff in future time steps
    DO JTIM=2,MERGE(DPARAM%NTDH_NEED,NTDH,USE_NTDH_NEED) ! update and move array of states within the routing convolution
      FUTURE(JTIM-1) = FUTURE(JTIM) + MROUTE%Q_INSTNT * DPARAM%FRAC_FUTURE(JTIM)
    END DO
    FUTURE(MTDH) = 0._wp  ! last element (just in case) - the rest are never accessed (treated as 0)
   
   CASE (iopt_no_routing)           ! no routing
    MROUTE%Q_ROUTED   = MROUTE%Q_INSTNT
   
  CASE DEFAULT                      ! check for errors
    message=trim(message)//"SMODL%iQ_TDH must be either iopt_rout_gamma or iopt_no_routing"
    ierr=10; return
  
  END SELECT
  ! ---------------------------------------------------------------------------------------
  END SUBROUTINE Q_OVERLAND

end module Q_OVERLAND_module
