module alloc_domain_module

  USE nrtype
  USE info_types, only: fuse_info
  USE data_types, only: domain_data

  implicit none
  private
  public :: allocate_domain_data
  public :: set_legacy_arrays

CONTAINS

  subroutine allocate_domain_data(info, domain, ierr, message)

  implicit none

  type(fuse_info),   intent(in)    :: info
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

  ! allocate ancillary forcing
  allocate(domain%ancil(nx,ny), stat=ierr)
  if(ierr/=0)then; message=trim(message)//"cannot allocate ancil"; return; endif

  ! allocate forcing window
  allocate(domain%force(nx,ny,nt), stat=ierr)
  if(ierr/=0)then; message=trim(message)//"cannot allocate force"; return; endif

  ! allocate state window
  allocate(domain%state(nx,ny,nt+1), stat=ierr)
  if(ierr/=0)then; message=trim(message)//"cannot allocate state"; return; endif

  ! allocate flux window
  allocate(domain%flux(nx,ny,nt), stat=ierr)
  if(ierr/=0)then; message=trim(message)//"cannot allocate flux"; return; endif

  ! allocate basin averages
  allocate(domain%aForce(nt), domain%aRoute(nt), stat=ierr)
  if(ierr/=0)then; message=trim(message)//"cannot allocate aForce/aRoute"; return; endif

  ! allocate routing if needed
  allocate(domain%route(nx,ny,nt), stat=ierr)
  if(ierr/=0)then; message=trim(message)//"cannot allocate route"; return; endif

  ! allocate bands
  allocate(domain%bands(nx,ny,nb,nt+1), stat=ierr)
  if(ierr/=0)then; message=trim(message)//"cannot allocate bands"; return; endif

  end subroutine allocate_domain_data

  ! -------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------

  ! ----- copy arrays in the domain structure to legacy arrays ---------------------

  subroutine set_legacy_arrays(info, domain, ierr, message)

  ! legacy modules
  use multiforce, only: nSpat1, nSpat2, numtim_sub
  use multiForce, only: gForce_3d, ancilF, aValid
  use multiState, only: gState_3d
  use multiRoute, only: AROUTE_3d
  use multiBands, only: MBANDS_VAR_4d, N_BANDS
  implicit none

  type(fuse_info),   intent(in)    :: info
  type(domain_data), intent(inout) :: domain
  integer(i4b),      intent(out)   :: ierr
  character(*),      intent(out)   :: message

  ierr = 0
  message = 'set_legacy_arrays/'

  ! ensure the spatial dimensions match what is in domain%info
  ! NOTE: dimension variables stored in legacy data modules
  nSpat1        = info%space%nx_local  ! NOTE: local to rank (MPI parallelization)
  nSpat2        = info%space%ny_local
  numtim_sub    = info%time%nt_window
  n_bands       = info%snow%n_bands

  ! allocate space for the forcing grid and states
  allocate(ancilF(nspat1,nspat2), stat=ierr)
  if(ierr/=0)then
    message=trim(message)//'unable to allocate space for 2-d ancillary grid'
    ierr=20; return
  endif
  
  ! allocate space for the forcing grid and states with a time dimension - only for subperiod
  allocate(aValid(   nspat1,nspat2,numtim_sub),   &
           gForce_3d(nspat1,nspat2,numtim_sub),   &
           gState_3d(nspat1,nspat2,numtim_sub+1), &
           AROUTE_3d(nspat1,nspat2,numtim_sub),   stat=ierr)
  if(ierr/=0)then
    message=trim(message)//'unable to allocate space for 3d structure'
    ierr=20; return
  endif
  
  ! allocate space for elevation bands
  allocate(MBANDS_VAR_4d(nspat1,nspat2,N_BANDS,numtim_sub+1), stat=ierr)
  if(ierr/=0)then
    message=trim(message)//'unable to allocate space for elevation bands'
    ierr=20; return
  endif

  ! copy arrays in the domain structure to legacy arrays
  aValid        = domain%valid   ! validity mask
  ancilF        = domain%ancil   ! ancillary forcing
  gForce_3d     = domain%force   ! forcing window
  gState_3d     = domain%state   ! state window
  AROUTE_3d     = domain%route   ! routing window
  MBANDS_VAR_4d = domain%bands   ! elevation band window

  end subroutine set_legacy_arrays

end module alloc_domain_module
