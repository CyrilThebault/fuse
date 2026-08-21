module time_windows_module

  use nrtype
  use info_types, only: fuse_info
  use time_utils, only: date_extractor, juldayss

  implicit none

  private
  public :: get_time_windows
  public :: export_time_to_multiforce

  contains

  subroutine get_time_windows(ncid, info, ierr, message)

    integer(i4b),      intent(in)    :: ncid
    type(fuse_info),   intent(inout) :: info
    integer(i4b),      intent(out)   :: ierr
    character(*),      intent(out)   :: message

    integer(i4b) :: nt
    character(len=1024) :: units_local
    integer(i4b) :: ios
    character(len=1024) :: cmessage

    ierr=0; message="get_time_windows/"

    ! ----- read forcing time axis ------------------------------------------------------

    call read_time_axis(ncid,                         &
                        info%time%time_steps,         &
                        info%time%time_bounds,        &
                        units_local, nt, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    info%time%nt_global = nt
    info%time%units     = trim(units_local)

    ! ----- build julian-day axis -------------------------------------------------------

    call build_julian_axis(info%time%time_steps,      &
                           info%time%time_bounds,     &
                           trim(units_local),         &
                           info%time%jdate_ref,       &
                           info%time%jdate,           &
                           info%time%deltim_days,     &
                           ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! ----- compute indices for sim/eval windows ----------------------------------------

    ! simulation indices
    call map_dates_to_indices(info%time%jdate, info%config%date_start_sim, info%config%date_end_sim, &
                              info%time%sim_beg, info%time%sim_end, ierr, cmessage)
    if (ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! evaluation indices
    call map_dates_to_indices(info%time%jdate, info%config%date_start_eval, info%config%date_end_eval, &
                              info%time%eval_beg, info%time%eval_end, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! ----- validate window consistency -------------------------------------------------

    call validate_windows(info%time, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! ----- derive simulation length ----------------------------------------------------

    info%time%nt_sim = info%time%sim_end - info%time%sim_beg + 1

    ! ----- configure sub-period windowing ----------------------------------------------

    ! convert sub-period string to integer
    read(info%config%numtim_sub_str,*,iostat=ios) info%time%nt_window
    if(ios/=0) then
      ierr=1; message=trim(message)//"cannot parse numtim_sub_str"; return
    endif

    ! handle cases where sub-periods are undefined
    if(info%time%nt_window == -9999) then
      info%time%use_subperiods = .false.
      info%time%nt_window      = info%time%nt_sim
    else
      info%time%use_subperiods = .true.
      ! keep nt_window as user-chosen chunk size
    endif

    ! ----- export info to legacy data structures ---------------------------------------

    ! export info%time -> multiforce to keep legacy code working
    call export_time_to_multiforce(info)

  end subroutine get_time_windows

  ! -------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------

  ! ----- backwards compatibility: export to multiforce globals -------------------------

  ! - New code stores all time-window metadata in info%time (source of truth).
  ! - Legacy routines still read multiforce globals (sim_beg, sim_end, numtim_sub, ...).

  subroutine export_time_to_multiforce(info)
    use multiforce, only: time_steps, timeUnits
    use multiforce, only: sim_beg, sim_end, eval_beg, eval_end, numtim_sim, numtim_sub, &
                          SUB_PERIODS_FLAG, istart, deltim
    implicit none
    type(fuse_info), intent(in) :: info

    time_steps = info%time%time_steps
    timeUnits  = info%time%units

    sim_beg    = info%time%sim_beg
    sim_end    = info%time%sim_end
    eval_beg   = info%time%eval_beg
    eval_end   = info%time%eval_end

    numtim_sim = info%time%nt_sim
    numtim_sub = info%time%nt_window
    SUB_PERIODS_FLAG = info%time%use_subperiods

    istart = sim_beg

    deltim = info%time%deltim_days
  end subroutine

  ! -------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------
  ! ----- helper routines ---------------------------------------------------------------
  ! -------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------

  ! ----- helper: read time axis from NetCDF --------------------------------------------

  subroutine read_time_axis(ncid, time_steps, time_bounds, units, nt, ierr, message)

    use netcdf

    implicit none

    integer(i4b),          intent(in)  :: ncid
    real(wp), allocatable, intent(out) :: time_steps(:)
    real(wp), allocatable, intent(out) :: time_bounds(:,:)
    character(len=*),      intent(out) :: units
    integer(i4b),          intent(out) :: nt
    integer(i4b),          intent(out) :: ierr
    character(*),          intent(out) :: message

    integer(i4b)       :: varid, bnd_varid
    integer(i4b)       :: dimids(1)
    character(len=256) :: bounds_name


    ierr=0; message="read_time_axis/"

    ! time dimension info

    ierr = nf90_inq_varid(ncid, "time", varid)
    if(ierr/=nf90_noerr) then
      message=trim(message)//"cannot find time variable"; return
    endif

    ierr = nf90_inquire_variable(ncid, varid, dimids=dimids)
    if(ierr/=nf90_noerr) then
      message=trim(message)//trim(nf90_strerror(ierr)); return
    endif

    ierr = nf90_inquire_dimension(ncid, dimids(1), len=nt)
    if(ierr/=nf90_noerr) then
      message=trim(message)//trim(nf90_strerror(ierr)); return
    endif

    ierr = nf90_get_att(ncid, varid, "units", units)
    if(ierr/=nf90_noerr) then
      message=trim(message)//"cannot read time units attribute"; return
    endif

    ! allocate space

    allocate(time_steps(nt), stat=ierr)
    if(ierr/=0) then
      message=trim(message)//"allocate(time_steps) failed"; return
    endif

    allocate(time_bounds(2,nt), stat=ierr)
    if(ierr/=0) then
      message=trim(message)//"allocate(time_bounds) failed"; return
    endif
   
    ! time bounds

    ierr = nf90_get_att(ncid, varid, "bounds", bounds_name)
    if(ierr/=nf90_noerr) then
      message=trim(message)//"time variable must define a bounds attribute"; return
    endif
   
    ierr = nf90_inq_varid(ncid, trim(bounds_name), bnd_varid)
    if(ierr/=nf90_noerr) then
      message=trim(message)//"cannot find time bounds variable"; return
    endif
   
    ierr = nf90_get_var(ncid, bnd_varid, time_bounds)
    if(ierr/=nf90_noerr) then
      message=trim(message)//"cannot read time bounds"; return
    endif

    ! FUSE convention: time is the midpoint of the interval
    time_steps = 0.5_wp * (time_bounds(1,:) + time_bounds(2,:))

  end subroutine read_time_axis

  ! -------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------

  ! ----- helper: build julian axis -----------------------------------------------------

  subroutine build_julian_axis(time_steps, time_bounds, units, jref, jdate, deltim_days, ierr, message)

    real(wp), intent(in) :: time_steps(:)
    real(wp), intent(in) :: time_bounds(:,:)
    character(len=*), intent(in) :: units
    real(wp), intent(out) :: jref
    real(wp), allocatable, intent(out) :: jdate(:)
    real(wp), intent(out)     :: deltim_days
    integer(i4b), intent(out) :: ierr
    character(*), intent(out) :: message

    integer(i4b) :: iy,im,id,ih
    character(len=1024) :: cmessage
    real(wp) :: scale_to_days

    real(wp) :: dt_current
    real(wp) :: tolerance
    integer(i4b) :: i
    logical(lgt), parameter :: do_timeCheck = .true.

    ierr=0; message="build_julian_axis/"

    ! extract reference date from the units string
    call date_extractor(trim(units), iy, im, id, ih)
    call juldayss(iy,im,id,ih, jref, ierr, cmessage)
    if(ierr/=0) then;  message=trim(message)//trim(cmessage); return; endif

    ! determine scaling factor to convert time_steps into days
    scale_to_days = time_units_to_days(units, ierr, cmessage)
    if(ierr/=0) then;  message=trim(message)//trim(cmessage); return; endif

    ! build julian axis (units of days)
    allocate(jdate(size(time_steps)), stat=ierr)
    if(ierr/=0) then; message=trim(message)//"allocate(jdate) failed"; return; endif
    jdate = jref + time_steps * scale_to_days

    ! define the forcing timestep length in days
    if (size(jdate) < 2) then
      ierr = 1
      message = trim(message)//"at least two forcing time steps are required"
      return
    endif

    deltim_days = (time_bounds(2,1) - time_bounds(1,1)) * scale_to_days

    if (deltim_days <= 0._wp) then
      ierr = 1
      message = trim(message)//"forcing time step must be positive"
      return
    endif

    ! Allow for floating-point round-off in converted time coordinates
    tolerance = epsilon(deltim_days) * 100._wp

    ! Verify that the forcing time axis is increasing and regularly spaced.
    if (do_timeCheck) then

      do i = 1, size(time_steps)

        dt_current = (time_bounds(2,i) - time_bounds(1,i)) * scale_to_days

        if (dt_current <= 0._wp) then
          ierr = 1
          message = trim(message)//"forcing time axis must be strictly increasing"
          return
        endif

        if (abs(dt_current - deltim_days) > tolerance) then
          ierr = 1
          message = trim(message)//"forcing time steps are not equally spaced"
          return
        endif

      enddo

    endif

  end subroutine build_julian_axis

  ! -------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------

  ! ----- helper: determine scaling factor to convert time_steps into days --------------

  real(wp) function time_units_to_days(units, ierr, message)
    implicit none
    character(len=*), intent(in) :: units
    integer(i4b), intent(out) :: ierr
    character(*), intent(out) :: message

    character(len=:), allocatable :: u
    integer(i4b) :: p

    ierr=0; message="time_units_to_days/"

    ! lower-case copy (simple approach)
    u = tolower_str( trim(adjustl(units)) )

    ! Look at the first token before a space
    p = index(u, " ")
    if(p <= 1) then
      ierr=1; message=trim(message)//"cannot parse units string: "//trim(units)
      time_units_to_days = 0._wp
      return
    endif

    select case (trim(u(1:p-1)))
      case ("days", "day")
        time_units_to_days = 1._wp
      case ("hours", "hour")
        time_units_to_days = 1._wp / 24._wp
      case ("minutes", "minute", "mins", "min")
        time_units_to_days = 1._wp / 1440._wp
      case ("seconds", "second", "secs", "sec")
        time_units_to_days = 1._wp / 86400._wp
      case default
        ierr=1
        message=trim(message)//"unsupported time unit: "//trim(u(1:p-1))
        time_units_to_days = 0._wp
    end select

  end function time_units_to_days

  pure function tolower_str(s) result(out)
    character(len=*), intent(in) :: s
    character(len=len(s)) :: out
    integer :: i
    do i=1,len(s)
      select case(s(i:i))
        case("A":"Z"); out(i:i) = achar(iachar(s(i:i)) + 32)
        case default;  out(i:i) = s(i:i)
      end select
    end do
  end function tolower_str

  ! -------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------

  ! ----- helper: map start/end date strings to indices ---------------------------------

  subroutine map_dates_to_indices(jdate, date_start, date_end, i_beg, i_end, ierr, message)

    real(wp), intent(in) :: jdate(:)
    character(len=*), intent(in) :: date_start, date_end
    integer(i4b), intent(out) :: i_beg, i_end
    integer(i4b), intent(out) :: ierr
    character(*), intent(out) :: message

    integer(i4b) :: iy,im,id,ih
    real(wp) :: j_start, j_end
    character(len=1024) :: cmessage

    ierr=0; message="map_dates_to_indices/"

    ! start date
    call date_extractor(trim(date_start), iy,im,id,ih)
    call juldayss(iy,im,id,ih, j_start, ierr, cmessage)
    if(ierr/=0) then;  message=trim(message)//trim(cmessage); return; endif

    ! end date
    call date_extractor(trim(date_end), iy,im,id,ih)
    call juldayss(iy,im,id,ih, j_end, ierr, cmessage)
    if(ierr/=0) then;  message=trim(message)//trim(cmessage); return; endif

    ! validate

    if(j_start > j_end) then
      ierr=1; message=trim(message)//"start date > end date"; return
    endif

    if(j_start < minval(jdate) .or. j_end > maxval(jdate)) then
      ierr=1; message=trim(message)//"requested window outside forcing range"; return
    endif

    ! get indices in jdate vector
    i_beg = minloc(abs(jdate - j_start), 1)
    i_end = minloc(abs(jdate - j_end  ), 1)

  end subroutine map_dates_to_indices


  ! -------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------

  ! ----- helper: validate sim/eval logic -----------------------------------------------

  subroutine validate_windows(ti, ierr, message)

    use info_types, only: time_info
    type(time_info), intent(in) :: ti
    integer(i4b), intent(out) :: ierr
    character(*), intent(out) :: message

    ierr=0; message="validate_windows/"

    if(ti%eval_beg < ti%sim_beg) then
      ierr=1; message=trim(message)//"eval start < sim start"; return
    endif
    if(ti%eval_end > ti%sim_end) then
      ierr=1; message=trim(message)//"eval end > sim end"; return
    endif

  end subroutine validate_windows

  ! -------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------

end module time_windows_module
