module physics_orig_module

  use nrtype
  use domain_types, only: domain_data
  use work_types,   only: fuse_work

  implicit none

  private
  public :: physics_orig

contains

  subroutine physics_orig(work,                     &
                          dt_sub, dt_full,          &
                          sub_idx, iSpat1, iSpat2,  &
                          ierr, message)

  ! shared data
  use multiforce,   only: DELTIM, MFORCE
  use multistate,   only: FSTATE, MSTATE
  use multiroute,   only: MROUTE
  use multibands,   only: MBANDS, Z_FORCING
  use multi_flux,   only: W_FLUX, M_FLUX

  ! original solver interface
  use interfaceb, only: ode_int, fuse_solve

  ! model options
  use model_defn,   only: SMODL
  use model_defnames

  implicit none

  type(fuse_work)       , intent(inout) :: work           ! work structures that depend on npar/nState
  integer(i4b)          , intent(in)    :: sub_idx, iSpat1, iSpat2
  real(wp)              , intent(inout) :: dt_sub, dt_full
  
  integer(i4b)          , intent(out)   :: ierr
  character(len=*)      , intent(out)   :: message

  ! locals
  character(len=1024)       :: cmessage

  ierr = 0
  message = "physics_orig/"

  ! copy data into legacy structures
  MFORCE = work%step%force
  W_FLUX = work%step%flux
  FSTATE = work%step%state0
  
  M_FLUX = W_FLUX
  MSTATE = FSTATE

  ! -------------------------
  ! snow module
  ! -------------------------
  select case(SMODL%iSNOWM)

    case (iopt_temp_index)

      ! copy data into legacy structures
      Z_FORCING      = work%snow%z_forcing
      MBANDS(:)%info = work%snow%sbands(:)%info
      MBANDS(:)%var  = work%snow%sbands(:)%var%bands_var

      ! run snow model
      call UPDATE_SWE(DELTIM)

      ! copy data back into canonical strucures
      work%snow%sbands(:)%var%bands_var = MBANDS(:)%var

    case (iopt_no_snowmod)
      call QRAINERROR()

    case default
      message=trim(message)//'unknown SMODL%iSNOWM option'
      ierr=10; return

  end select 

  ! -------------------------
  ! interception
  ! -------------------------
  
  M_FLUX%PIN0 = M_FLUX%EFF_PPT
  call UPDATE_INTERCEPTION(DELTIM, ierr, cmessage)
  if (ierr /= 0) then; message = trim(message)//trim(cmessage); return; end if

  ! -------------------------
  ! soil physics
  ! -------------------------
  
  call ODE_INT(FUSE_SOLVE, work%num%x0, work%num%x1, dt_sub, dt_full, ierr, cmessage)
  if (ierr /= 0) then; message=trim(message)//trim(cmessage); return; end if

  ! save fluxes in canonical structures
  work%step%flux  = W_FLUX

  end subroutine physics_orig

end module physics_orig_module
