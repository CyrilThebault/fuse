module setup_domain_module

  USE nrtype
  USE info_types, only: cli_options
  USE info_types, only: fuse_info
  USE data_types, only: domain_data
  USE fuse_globaldata, only: isPrint

  implicit none

  private
  public :: setup_domain

contains

  subroutine setup_domain(opts, info, domain, ierr, message)

  ! access subroutines
  use netcdf, only: nf90_open, nf90_nowrite, nf90_strerror      ! NetCDF functions
  USE fuse_fileManager,      only: read_fuse_control_file       ! sets directories and filenames

  USE domain_dims_module,    only: get_domain_dims              ! get nx, ny, nt, and nbands
  USE domain_decomp_module,  only: get_domain_decomp_indices    ! get MPI domain decomposition indices
  
  USE time_windows_module,   only: get_time_windows             ! get info on the rolling time windows
  USE time_windows_module,   only: export_time_to_multiforce    ! populate legacy multiforce modules

  USE get_gforce_module,     only: get_forcing_varids           ! get name/varid table for forcing variables
  USE get_gforce_module,     only: read_latlon_2d               ! read lat/lon
  USE read_elevbands_module, only: read_elevbands               ! read elevation bands

  USE alloc_domain_module,   only: allocate_domain_data         ! allocate space for data arrays in the domain structure
  USE alloc_domain_module,   only: set_legacy_arrays            ! copy arrays in the domain%data structure to legacy arrays 

  USE init_mizuRoute_topo,   only: init_mizuroute_topology      ! initialize mizuRoute network topology

  implicit none
  
  ! input
  type(cli_options)   , intent(in)                  :: opts            ! command line interface options
  type(fuse_info)     , intent(inout)               :: info            ! domain info
  type(domain_data)   , intent(inout)               :: domain          ! domain data
  
  ! output
  integer(i4b)        , intent(out)                 :: ierr            ! error code
  character(len=1024) , intent(out)                 :: message         ! error message
  
  ! ----- internal -----------------------------------------------------------------------
  CHARACTER(LEN=1024)                               :: CMESSAGE        ! error message
  ! ---------------------------------------------------------------------------------------
  ierr=0; message='setup_domain/'

  ! ----- set paths and file names --------------------------------------------------------
  
  ! read fuse control file (set paths/filenames etc.)
  call read_fuse_control_file(trim(opts%control_file), opts, info, ierr, cmessage) 
  if (ierr/=0)then; message=trim(message)//trim(cmessage); ierr=20; return; endif

  ! ----- initialize the river-network topology ------------------------------------------

  call init_mizuroute_topology(info, domain, ierr, cmessage)
  if (ierr/=0)then; message=trim(message)//trim(cmessage); ierr=20; return; endif

  ! ----- read domain metadata ------------------------------------------------------------
  
  ! populate domain structure with dimension lengths
  !   -- nx_global, ny_global, nt_global, n_bands
  call get_domain_dims(info, ierr, cmessage)
  if (ierr/=0)then; message=trim(message)//trim(cmessage); ierr=20; return; endif

  ! get indices for MPI decomposition of the spatial domain: y_start_global, ny_local 
  ! NOTE: These indices will be used later to read different subsets of forcing data for different ranks
  call get_domain_decomp_indices(info)

  ! ----- read grid info and define indices for MPI domain decomposition ------------------

  ! open NetCDF forcing file
  ierr = nf90_open(trim(info%files%fname_netcdf_forc), nf90_nowrite, info%files%ncid_forc)
  if (ierr/=0)then; message=trim(message)//' nf90_open failed: '//trim(nf90_strerror(ierr)); return; endif
  if(isPrint) print *, 'Open forcing file:', trim(info%files%fname_netcdf_forc)
  if(isPrint) PRINT *, 'NCID_FORC is', info%files%ncid_forc
 
  ! ----- Compute time indices for sim/eval windows and subperiod chunk size --------------
  !
  ! Reads the forcing-file NetCDF time coordinate (and units), and:
  !   - builds a Julian-day time axis and timestep in days
  !   - maps the user-specified simulation/evaluation date ranges into index windows
  !     (and optional subperiod chunks) stored in info%time
  call get_time_windows(info%files%ncid_forc, info, ierr, cmessage)
  if (ierr/=0)then; message=trim(message)//trim(cmessage); ierr=20; return; endif

  ! ----- Allocate space for domain data --------------------------------------------------

  ! allocate space for the arrays in the domain%data structure
  call allocate_domain_data(info, domain, ierr, cmessage)
  if (ierr/=0)then; message=trim(message)//trim(cmessage); ierr=20; return; endif

  ! ----- Read lat/lon, elevation band arrays, and forcing var ids  -----------------------

  ! read lat/lon and store in the domain%data%coords structure
  call read_latlon_2d(info%files%ncid_forc, info, domain%coords, ierr, message)
  if (ierr/=0)then; message=trim(message)//trim(cmessage); ierr=20; return; endif

  ! read elevation bands information and store in the domain%data structure
  call read_elevbands(info, domain, ierr, cmessage)
  if (ierr/=0)then; message=trim(message)//trim(cmessage); ierr=20; return; endif

  call get_forcing_varids(info%files%ncid_forc, info, ierr, cmessage)
  if (ierr/=0)then; message=trim(message)//trim(cmessage); ierr=20; return; endif

  ! ----- Routines that use the old structures --------------------------------------------

  ! copy arrays in the domain structure to legacy arrays
  call set_legacy_arrays(info, domain, ierr, cmessage)
  if (ierr/=0)then; message=trim(message)//trim(cmessage); ierr=20; return; endif

  ! defines method/parameters used for numerical solution based on numerix file
  ! NOTE: This routine supports the legacy FUSE v1 numerics experiments
  CALL GETNUMERIX(IERR,CMESSAGE)
  if (ierr/=0)then; message=trim(message)//trim(cmessage); ierr=20; return; endif

  print*, 'end of setup_domain'

  end subroutine setup_domain

  ! -------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------

end module setup_domain_module
