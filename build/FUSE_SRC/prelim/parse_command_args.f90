module parse_command_args_MODULE

  USE nrtype
  USE info_types, only: cli_options

  implicit none

  private
  public :: parse_command_args

contains

  subroutine parse_command_args(opts, err, message)
  implicit none
  ! dummies
  type(cli_options)   , intent(out)   :: opts          ! command line interface options
  integer(i4b)        , intent(out)   :: err           ! error code
  character(len=1024) , intent(out)   :: message       ! error message
  ! internal
  integer(i4b)                        :: i             ! index of command line argument
  character(len=:)    , allocatable   :: a, v          ! command line arguments
  character(len=:)    , allocatable   :: cIndex        ! character index
  character(len=:)    , allocatable   :: kv, pname, pval_str  ! parameter strings
  real(sp)                            :: pval          ! parameter value
  integer(i4b)                        :: nArg          ! number of command line arguments
  character(len=:)    , allocatable   :: argname
  character(len=:)    , allocatable   :: cmessage
  ! initialize error control
  err=0; message='parse_command_args/'
  
  ! ----- parse command line arguments ------------------------------------------------------

  ! -----------------------------------------------------------------------------------------
  ! Parse command-line options
  !   -c, --control   <file>             required unless --help/--version
  !   -m, --runmode   <def|idx|opt|sce>  required unless --help/--version
  !   -d, --domid     <string>           required unless --help/--version
  !   -s, --sets      <file>             required for idx,opt
  !   -i, --index     <int>              required for idx
  !   -r, --restart   <y|m|d|e|never>    optional
  !   -t, --tag       <string>           optional
  !   -p, --param     NAME=VALUE         repeatable
  !   -v, --version                      print version and exit
  !   -h, --help                         print help and exit
  ! -----------------------------------------------------------------------------------------

  nArg = command_argument_count()

  i = 1
  do while (i <= narg)
    call get_arg(i,a)

    select case (trim(a))

    case ('-h','--help')
      opts%show_help = .true.
      i = i + 1

    case ('-v','--version')
      opts%show_version = .true.
      i = i + 1

     case ('-t','--tag')
      call require_next(i, narg, a, v, err, cmessage)
      opts%tag = trim(v)
      i = i + 2

    case ('-c','--control')
      call require_next(i, narg, a, v, err, cmessage)
      opts%control_file = trim(v)
      i = i + 2

    case ('-m','--runmode')
      call require_next(i, narg, a, v, err, cmessage)
      opts%runmode = to_lower(trim(v))
      i = i + 2

    case ('-d','--domid')
      call require_next(i, narg, a, v, err, cmessage)
      opts%domain_id = trim(v)
      i = i + 2

    case ('-p', '--param')
      call require_next(i, narg, a, kv, err, cmessage)
      i = i + 2

    case ('-s','--sets','--param-sets')
      call require_next(i, narg, a, v, err, cmessage)
      opts%sets_file = trim(v)
      i = i + 2

    case ('-i','--index')
      call require_next(i, narg, a, cIndex, err, cmessage)
      i = i + 2

    case ('-r','--restart')
      call require_next(i, narg, a, v, err, cmessage)
      opts%restart_freq = to_lower(trim(v))
      i = i + 2

    case default
      if (len_trim(a) > 0 .and. a(1:1) == '-') then
        err = 1
        cmessage = "unknown option: "//trim(a)//"; type 'fuse.exe --help' for usage"
      else
        err = 1
        cmessage = "unexpected positional argument: "//trim(a)//"; type 'fuse.exe --help' for usage"
      end if
    end select

    ! process error code
    if(err/=0)then
     message=trim(message)//trim(cmessage)
     err=20; return
    endif

    ! process parameters -- needs to be in the do loop since multiple parameters
    if(allocated(kv))then

      ! split name/value based on the equal sign
      call split_param_kv(trim(kv), pname, pval_str, err, cmessage)
      if(err /= 0)then; message=trim(message)//trim(cmessage); err=20; return; endif

      ! convert characters to real values
      call parse_real_sp(pval_str, pval, err, cmessage)
      if (err /= 0) then
        message=trim(message)//"invalid --param value for "//trim(pname)//": "//trim(cmessage)
        err=20; return
      end if

      ! add to structure in opts
      call push_param(opts%param_name, opts%param_value, pname, pval)
      print*, opts%param_name
      print*, opts%param_value
    
    endif  ! if processing parameters

  end do  ! looping through arguments

  ! Early exits
  if (opts%show_help) then
    call printCommandHelp()
    stop 0
  end if
  if (opts%show_version) then
    call printVersionInfo()
    stop 0
  end if

  ! Parse parameter index
  if(allocated(cIndex))then
   call parse_int(cIndex, opts%indx, err, cmessage)
   if(err/=0)then
     message=trim(message)//trim(cmessage)
     err=20; return
    endif
  endif

  ! Validate required args
  
  argname = ""
  
  call check_required(allocated(opts%domain_id),    "--domid",   argname)
  call check_required(allocated(opts%control_file), "--control", argname)
  call check_required(allocated(opts%runmode),      "--runmode", argname)
  
  if (len_trim(argname) > 0) then
    err = 1
    message = trim(message)//"Missing required arguments:"//trim(argname)// &
              new_line("a")//new_line("a")// &
              "Usage:"//new_line("a")// &
              "  fuse.exe -d DOMAIN -c CONTROL -m MODE [options]"// &
              new_line("a")//new_line("a")// &
              "Try 'fuse.exe --help' for detailed usage."
    return
  end if

  ! Check run mode
  if (.not. is_valid_mode(opts%runmode)) then
    err = 1; message = trim(message)//"invalid --runmode: "//trim(opts%runmode)//" (expect def|idx|opt|sce)"; return
  end if

  ! Mode-dependent requirements
  select case (trim(opts%runmode))
  case ('idx')
    if (.not. allocated(opts%sets_file)) then
      err = 1; message = trim(message)//"runmode idx requires --sets <file>"; return
    end if
    if (opts%indx < 0) then
      err = 1; message = trim(message)//"runmode idx requires --index <int>"; return
    end if
  case ('opt')
    if (.not. allocated(opts%sets_file)) then
      err = 1; message = trim(message)//"runmode opt requires --sets <file>"; return
    end if
  case ('def','sce')
    ! no extra requirements
  end select

  ! Validate frequencies if provided (optional)
  if (allocated(opts%restart_freq)) then
    if (.not. is_valid_restart(opts%restart_freq)) then
      err = 1; message = trim(message)//"invalid --restart: "//trim(opts%restart_freq)//" (expect y|m|d|e|never)"; return
    end if
  end if

  contains

    subroutine check_required(is_present, name, message)
      logical, intent(in) :: is_present
      character(*), intent(in) :: name
      character(:), allocatable, intent(inout) :: message
    
      if (.not. is_present) then
        message = trim(message)//new_line("a")//"  "//name
      end if
    end subroutine check_required

  end subroutine parse_command_args

  ! ----- list version ----------------------------------------------------------------------
  
  subroutine printVersionInfo()
    ! Assumes these are available, e.g. from:
    !   include "fuseversion.inc"
    ! somewhere in a used module (e.g., fuse_globaldata) OR add that include here.
    use fuse_globaldata, only: FUSE_VERSION, FUSE_BUILDTIME, FUSE_GITBRANCH, FUSE_GITHASH
    implicit none
    print '(A)', repeat('-', 70)
    print '(A)', 'FUSE'
    print '("  ",A12," : ",A)', 'Version',    trim(FUSE_VERSION)
    print '("  ",A12," : ",A)', 'Build time', trim(FUSE_BUILDTIME)
    print '("  ",A12," : ",A)', 'Git branch', trim(FUSE_GITBRANCH)
    print '("  ",A12," : ",A)', 'Git hash',   trim(FUSE_GITHASH)
    print '(A)', repeat('-', 70)
  end subroutine printVersionInfo

  ! ----- list command usage ----------------------------------------------------------------

  subroutine printCommandHelp()
    implicit none
    print "(A)", ""
    print "(A)", "FUSE: Framework for Understanding Structural Errors"
    print "(A)", ""
    print "(A)", "Usage:"
    print "(A)", "  fuse.exe -d DOMAIN -c CONTROL -m MODE [options]"
    print "(A)", ""
    print "(A)", "Required:"
    print "(A)", "  -d, --domid <string>        Domain identifier"
    print "(A)", "  -c, --control <file>        FUSE control file"
    print "(A)", "  -m, --runmode <mode>        def | idx | opt | sce"
    print "(A)", ""
    print "(A)", "Run modes:"
    print "(A)", "  def   Run with default parameters"
    print "(A)", "  idx   Run parameter-set INDEX from --sets"
    print "(A)", "  opt   Run best parameter set from --sets"
    print "(A)", "  sce   Calibrate parameters using SCE"
    print "(A)", ""
    print "(A)", "Conditional:"
    print "(A)", "  -s, --sets <file>           Parameter-set file (idx,opt)"
    print "(A)", "  -i, --index <integer>       Parameter-set index (idx)"
    print "(A)", ""
    print "(A)", "Optional:"
    print "(A)", "  -p, --param NAME=VALUE      Override parameter (repeatable)"
    print "(A)", "  -r, --restart <freq>        y | m | d | e | never"
    print "(A)", "  -t, --tag <string>          Append tag to output files"
    print "(A)", "  -v, --version               Print version and exit"
    print "(A)", "  -h, --help                  Print this help message"
    print "(A)", ""
    print "(A)", "Restart frequencies:"
    print "(A)", "  y=yearly  m=monthly  d=daily  e=every timestep"
    print "(A)", ""
    print "(A)", "Examples:"
    print "(A)", "  fuse.exe -d CAN_05BB001 -c control.txt -m def"
    print "(A)", "  fuse.exe -d CAN_05BB001 -c control.txt -m idx \\"
    print "(A)", "           -s params.nc -i 17"
    print "(A)", "  fuse.exe -d CAN_05BB001 -c control.txt -m opt \\"
    print "(A)", "           -s params.nc"
    print "(A)", "  fuse.exe -d CAN_05BB001 -c control.txt -m def \\"
    print "(A)", "           -p MAXWATR_1=1000 -p PERCRTE=0.25"
    print "(A)", ""
  end subroutine printCommandHelp

  ! -----------------------------------------------------------------------------------------
  ! Helpers
  ! -----------------------------------------------------------------------------------------

  subroutine get_arg(i, out)
    integer, intent(in) :: i
    character(len=:), allocatable, intent(out) :: out
    integer :: L
    call get_command_argument(i, length=L)
    allocate(character(len=L) :: out)
    call get_command_argument(i, out)
  end subroutine get_arg

  subroutine require_next(i, narg, opt, val, err, message)
    integer, intent(in) :: i, narg
    character(len=*), intent(in) :: opt
    character(len=:), allocatable, intent(out) :: val
    integer, intent(out) :: err
    character(len=:), allocatable, intent(out) :: message
    err = 0
    message = ""
    if (i+1 > narg) then
      err = 1
      message = "missing value after "//trim(opt)//"; type 'fuse.exe --help' for usage"
      return
    end if
    call get_arg(i+1, val)
  end subroutine require_next

  subroutine split_param_kv(kv, name, val, err, message)
    character(len=*), intent(in) :: kv
    character(len=:), allocatable, intent(out) :: name, val
    integer(i4b), intent(out) :: err
    character(len=:), allocatable, intent(out) :: message
    integer(i4b) :: p

    err = 0; message = ""
    p = index(kv, '=')
    if (p <= 1 .or. p >= len_trim(kv)) then
      err = 1
      message = "expected NAME=VALUE after --param, got: "//trim(kv)
      return
    end if

    name = adjustl(kv(1:p-1))
    val  = adjustl(kv(p+1:))

    if (len_trim(name) == 0 .or. len_trim(val) == 0) then
      err = 1
      message = "expected NAME=VALUE after --param, got: "//trim(kv)
      return
    end if
  end subroutine split_param_kv

  subroutine parse_real_sp(s, x, err, message)
    character(len=*), intent(in) :: s
    real(sp), intent(out) :: x
    integer, intent(out) :: err
    character(len=:), allocatable, intent(out) :: message
    integer(i4b) :: ios
    err = 0; message = ""
    read(s, *, iostat=ios) x
    if (ios /= 0) then
      err = 1
      message = "invalid real: "//trim(s)
    end if
  end subroutine parse_real_sp

  subroutine parse_int(s, x, err, message)
    character(len=*), intent(in) :: s
    integer, intent(out) :: x
    integer, intent(out) :: err
    character(len=:), allocatable, intent(out) :: message
    integer :: ios
    err = 0
    message = ""
    read(s, *, iostat=ios) x
    if (ios /= 0) then
      err = 1
      message = "invalid integer: "//trim(s)
    end if
  end subroutine parse_int

  pure function to_lower(s) result(t)
    character(len=*), intent(in) :: s
    character(len=len(s)) :: t
    integer :: k, c
    t = s
    do k = 1, len(s)
      c = iachar(t(k:k))
      if (c >= iachar('A') .and. c <= iachar('Z')) then
        t(k:k) = achar(c + (iachar('a') - iachar('A')))
      end if
    end do
  end function to_lower

  subroutine push_param(pnames, pvals, name, val)
    use nrtype
    implicit none
    character(len=:), allocatable, intent(inout) :: pnames(:)
    real(sp), allocatable, intent(inout)         :: pvals(:)
    character(len=*), intent(in)                 :: name
    real(sp), intent(in)                         :: val
   
    character(len=:), allocatable :: new_names(:)
    real(sp), allocatable         :: new_vals(:)
    integer :: n
   
    n = 0
    if (allocated(pvals)) n = size(pvals)
   
    allocate(character(len=len_trim(name)) :: new_names(n+1))
    allocate(new_vals(n+1))
   
    if (n > 0) then
      new_names(1:n) = pnames
      new_vals(1:n)  = pvals
    end if
   
    new_names(n+1) = trim(name)
    new_vals(n+1)  = val
   
    call move_alloc(new_names, pnames)
    call move_alloc(new_vals,  pvals)
  end subroutine push_param

  pure logical function is_valid_mode(m)
    character(len=*), intent(in) :: m
    is_valid_mode = (trim(m) == 'def' .or. trim(m) == 'idx' .or. trim(m) == 'opt' .or. trim(m) == 'sce')
  end function is_valid_mode

  pure logical function is_valid_restart(f)
    character(len=*), intent(in) :: f
    is_valid_restart = (trim(f) == 'y' .or. trim(f) == 'm' .or. trim(f) == 'd' .or. trim(f) == 'e' .or. trim(f) == 'never')
  end function is_valid_restart

end module parse_command_args_MODULE



