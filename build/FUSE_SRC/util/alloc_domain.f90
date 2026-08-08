module alloc_domain_module

  USE nrtype
  USE info_types, only: fuse_info
  USE domain_types, only: domain_data
 
  use fuse_globaldata, only: NVAR_FORC

  implicit none
  private
  public :: allocate_domain_data
  public :: set_legacy_arrays

CONTAINS

  subroutine allocate_domain_data(info, domain, ierr, message)

  implicit none

  type(fuse_info),   intent(inout) :: info
  type(domain_data), intent(inout) :: domain
  integer(i4b),      intent(out)   :: ierr
  character(*),      intent(out)   :: message

  integer(i4b) :: nx, ny, nt, nb

  ierr=0; message="allocate_domain_data/"

  ! define dimensions
  nx = info%space%nx_local  ! NOTE: local to rank (MPI parallelization)
  ny = info%space%ny_local
  nt = info%time%nt_window
  nb = info%snow%n_bands

  ! allocate validity mask
  allocate(domain%valid(nx,ny,nt), stat=ierr)
  if(ierr/=0)then; message=trim(message)//"cannot allocate valid"; return; endif

  ! allocate forcing window
  allocate(domain%force(nx,ny,nt), stat=ierr)
  if(ierr/=0)then; message=trim(message)//"cannot allocate force"; return; endif

  ! allocate state window
  allocate(domain%state(nx,ny,nt+1), stat=ierr)
  if(ierr/=0)then; message=trim(message)//"cannot allocate state"; return; endif

  ! allocate flux window
  allocate(domain%flux(nx,ny,nt), stat=ierr)
  if(ierr/=0)then; message=trim(message)//"cannot allocate flux"; return; endif

  ! allocate routing if needed
  allocate(domain%route(nx,ny,nt), stat=ierr)
  if(ierr/=0)then; message=trim(message)//"cannot allocate route"; return; endif

  ! allocate elevation grid
  allocate(domain%z_forcing(nx,ny), domain%elev_mask(nx,ny), stat=ierr)
  if(ierr/=0)then; message=trim(message)//"cannot allocate elev grid"; return; endif
  
  ! allocate elevation bands (info)
  allocate(domain%bands_info(nx,ny,nb), stat=ierr)
  if(ierr/=0)then; message=trim(message)//"cannot allocate elev bands (info)"; return; endif

  ! allocate elevation bands (var)
  allocate(domain%bands_var(nx,ny,nb,nt+1), stat=ierr)
  if(ierr/=0)then; message=trim(message)//"cannot allocate elev bands (var)"; return; endif

  ! allocate forcing lookup table
  allocate(info%files%forc%name(NVAR_FORC), info%files%forc%varid(NVAR_FORC), info%files%forc%multiplier(NVAR_FORC), stat=ierr)
  if(ierr/=0)then; message=trim(message)//"cannot allocate forcing lookup table"; return; endif
  info%files%forc%multiplier(:) = -1._wp

  end subroutine allocate_domain_data

  ! -------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------

  ! ----- copy arrays in the domain structure to legacy arrays ---------------------

  subroutine set_legacy_arrays(info, domain, ierr, message)

  ! legacy modules
  use multiforce, only: nSpat1, nSpat2, numtim_sub
  USE multiforce, only: startSpat2                         ! starting y index for data read
  USE multiforce, only: ncid_forc
  use multiforce, only: timeUnits
  use multiforce, only: nInput
  use multiForce, only: aValid
  use multiBands, only: N_BANDS, MBANDS
  implicit none

  type(fuse_info),   intent(in)    :: info
  type(domain_data), intent(inout) :: domain
  integer(i4b),      intent(out)   :: ierr
  character(*),      intent(out)   :: message

  integer(i4b) :: nx, ny, nt, nb

  ierr = 0
  message = 'set_legacy_arrays/'

  ! define dimensions
  nx = info%space%nx_local  ! NOTE: local to rank (MPI parallelization)
  ny = info%space%ny_local
  nt = info%time%nt_window
  nb = info%snow%n_bands

  ! set variables in multiforce
  nSpat1     = nx
  nSpat2     = ny
  numtim_sub = nt
  startSpat2 = info%space%y_start_global

  ncid_forc  = info%files%ncid_forc

  nInput     = info%config%nInput

  ! set bands
  N_BANDS = nb

  ! allocate validity mask
  allocate(aValid(nx,ny,nt), stat=ierr)
  if(ierr/=0)then; message=trim(message)//"cannot allocate valid"; return; endif

  ! allocate elevation bands
  allocate(MBANDS(nb), stat=ierr)
  if(ierr/=0)then; message=trim(message)//"cannot allocate elev bands"; return; endif

  ! copy arrays in the domain structure to legacy arrays
  aValid         = domain%valid      ! validity mask

  end subroutine set_legacy_arrays

end module alloc_domain_module
