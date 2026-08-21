module physics_diff_module

  use nrtype
  use domain_types, only: domain_data
  use work_types,   only: fuse_work

  implicit none

  private
  public :: physics_diff

contains

  subroutine physics_diff(work, dt_full,            &
                          sub_idx, iSpat1, iSpat2,  &
                          ierr, message)

  use implicit_solve_module,   only: implicit_solve
  use update_swe_diff_module,  only: update_swe_diff

  ! model options
  use model_defn,   only: SMODL, NSTATE
  use model_defnames

  implicit none

  type(fuse_work)       , intent(inout) :: work           ! work structures that depend on npar/nState
  real(wp)              , intent(in)    :: dt_full        ! time step length
  integer(i4b)          , intent(in)    :: sub_idx, iSpat1, iSpat2
  
  integer(i4b)          , intent(out)   :: ierr
  character(len=*)      , intent(out)   :: message

  ! locals
  character(len=1024)       :: cmessage

  ierr = 0
  message = "physics_diff/"

  ! -------------------------
  ! snow module
  ! -------------------------
  select case(SMODL%iSNOWM)

    case (iopt_temp_index)
      call UPDATE_SWE_DIFF(work, dt_full)

    case (iopt_no_snowmod)
      call QRAINERROR()

    case default
      message=trim(message)//'unknown SMODL%iSNOWM option'
      ierr=10; return

  end select 

  ! -------------------------
  ! interception
  ! -------------------------
  
  if (SMODL%iINTRC /= iopt_no_intrcep) then
    message = 'interception not yet implemented for differentiable mode'
    ierr = 1; return
  endif

  ! -------------------------
  ! soil physics
  ! -------------------------
  
  call implicit_solve(work, work%num%x0, work%num%x1, nState, ierr, cmessage)
  if (ierr /= 0) then; message=trim(message)//trim(cmessage); return; end if

  end subroutine physics_diff

end module physics_diff_module
