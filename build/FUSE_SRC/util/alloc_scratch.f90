module alloc_scratch_module


  USE nrtype
  use info_types, only: fuse_info
  use work_types, only: fuse_work

  implicit none
  private
  public :: init_fuse_work

CONTAINS

  subroutine init_fuse_work(info, work, ierr, message)

    use globaldata, only: NPAR_SNOW
    implicit none

    type(fuse_info),   intent(in)    :: info
    type(fuse_work),   intent(inout) :: work
    integer(i4b),      intent(out)   :: ierr
    character(*),      intent(out)   :: message

    integer(i4b) :: ib
    integer(i4b) :: nBands, nState, nPar
    logical(lgt) :: redo

    ierr=0; message="init_fuse_work/"

    nBands = info%snow%n_bands
    nState = info%config%nState
    nPar   = info%config%nParam

    ! check if there is a need to reallocate
    redo = needs_realloc_work(work, nBands, nState, nPar, NPAR_SNOW)
    if (.not. redo) return

    ! if need to reallocate, then need to free up space
    if (work%is_initialized) then
     call free_fuse_work(work, ierr, message)
     if(ierr/=0) return
    endif

    ! ---- allocate state vectors ----
    allocate(work%num%x0(nState), &
             work%num%x1(nState), stat=ierr)
    if(ierr/=0) then
      message=trim(message)//"cannot allocate state vectors"
      return
    endif

    ! optional debug scratch
    ! allocate(work%dSdt(nState), work%J(nState,nState), stat=ierr)

    ! ---- allocate differentiable parent derivatives ----
    allocate(work%adj%df_dS(nState), &
             work%adj%df_dPar(nPar), &
             work%adj%dL_dPar(nPar), stat=ierr)
    if(ierr/=0) then
      message=trim(message)//"cannot allocate derivatives"
      return
    endif

    ! ---- allocate elevation band containers ----
    allocate(work%snow%sbands(nBands), stat=ierr)
    if(ierr/=0) then
      message=trim(message)//"cannot allocate sbands"
      return
    endif

    ! ---- allocate per-band parameter derivative vectors ----
    do ib=1,nBands
      allocate(work%snow%sbands(ib)%var%dSWE_dParam(nPar_snow), stat=ierr)
      if(ierr/=0) then
        message=trim(message)//"cannot allocate dSWE_dParam for band"
        return
      endif
      work%snow%sbands(ib)%var%dSWE_dParam(:) = 0._sp
    enddo

    ! ---- initialize the band snow vars once ----
    work%snow%sbands(:)%var%SWE         = 0._sp
    work%snow%sbands(:)%var%SNOWACCMLTN = 0._sp
    work%snow%sbands(:)%var%SNOWMELT    = 0._sp
    work%snow%sbands(:)%var%DSWE_DT     = 0._sp

    work%is_initialized = .true.

  end subroutine init_fuse_work

  ! -------------------------------------------------------------------------------------

  subroutine free_fuse_work(work, ierr, message)

    implicit none
    type(fuse_work), intent(inout) :: work
    integer(i4b),    intent(out)   :: ierr
    character(*),    intent(out)   :: message

    integer(i4b) :: ib, istat

    ierr    = 0
    message = "free_fuse_work/"

    ! ---- state vectors ----
    if(allocated(work%num%x0)) then
      deallocate(work%num%x0, stat=istat)
      call note_fail("num%x0", istat)
    endif

    if(allocated(work%num%x1)) then
      deallocate(work%num%x1, stat=istat)
      call note_fail("num%x1", istat)
    endif

    ! ---- derivative arrays ----
    if (allocated(work%adj%df_dS)) then
      deallocate(work%adj%df_dS, stat=istat)
      call note_fail("adj%df_dS", istat)
    endif

    if (allocated(work%adj%df_dPar)) then
      deallocate(work%adj%df_dPar, stat=istat)
      call note_fail("adj%df_dPar", istat)
    endif

    if (allocated(work%adj%dL_dPar)) then
      deallocate(work%adj%dL_dPar, stat=istat)
      call note_fail("adj%dL_dPar", istat)
    endif

    ! ---- elevation band structures ----
    if (allocated(work%snow%sbands)) then

      do ib = 1, size(work%snow%sbands)
        if (allocated(work%snow%sbands(ib)%var%dSWE_dParam)) then
          deallocate(work%snow%sbands(ib)%var%dSWE_dParam, stat=istat)
          call note_fail("sbands%var%dSWE_dParam", istat)
        endif
      enddo

      deallocate(work%snow%sbands, stat=istat)
      call note_fail("snow%sbands", istat)

    endif

    work%is_initialized = .false.

    contains

      subroutine note_fail(where, istat)
        character(*), intent(in) :: where
        integer(i4b), intent(in) :: istat

        if (istat /= 0) then
          ! preserve the first nonzero stat as ierr
          if (ierr == 0) ierr = istat

          ! append context (do not overwrite)
          message = trim(message)//" dealloc_fail("//trim(where)//")"
        endif
      end subroutine note_fail

  end subroutine free_fuse_work

  ! -------------------------------------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------------------------------------

  ! Private: decide if we need to free+reallocate work arrays
  logical(lgt) function needs_realloc_work(work, nBands, nState, nPar, nParSnow) result(redo)

    type(fuse_work), intent(in) :: work
    integer(i4b),    intent(in) :: nBands, nState, nPar, nParSnow

    integer(i4b) :: ib

    redo = .false.

    ! Not initialized => must allocate
    if (.not. work%is_initialized) then
      redo = .true.
      return
    endif

    ! Must be allocated if we claim initialized
    if (.not. allocated(work%adj%df_dS))   then; redo=.true.; return; endif
    if (.not. allocated(work%adj%df_dPar)) then; redo=.true.; return; endif
    if (.not. allocated(work%adj%dL_dPar)) then; redo=.true.; return; endif
    if (.not. allocated(work%snow%sbands)) then; redo=.true.; return; endif

    ! Size checks
    if (size(work%adj%df_dS)    /= nState) then; redo=.true.; return; endif
    if (size(work%adj%df_dPar) /= nPar)   then; redo=.true.; return; endif
    if (size(work%adj%dL_dPar) /= nPar)   then; redo=.true.; return; endif
    if (size(work%snow%sbands) /= nBands) then; redo=.true.; return; endif

    ! Per-band arrays
    do ib = 1, nBands
      if (.not. allocated(work%snow%sbands(ib)%var%dSWE_dParam)) then
        redo = .true.; return
      endif
      if (size(work%snow%sbands(ib)%var%dSWE_dParam) /= nParSnow) then
        redo = .true.; return
      endif
    enddo

  end function needs_realloc_work


end module alloc_scratch_module
