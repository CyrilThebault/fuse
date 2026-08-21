module setup_domain_module

  USE nrtype
  USE info_types, only: cli_options
  USE info_types, only: fuse_info
  USE domain_types, only: domain_data
  USE fuse_globaldata, only: isPrint

  implicit none

  private
  public :: setup_domain

contains

  subroutine setup_domain(info, domain, ierr, message)

  ! access subroutines
  use netcdf, only: nf90_open, nf90_nowrite, nf90_strerror      ! NetCDF functions
  USE fuse_fileManager,      only: read_fuse_control_file       ! sets directories and filenames

  USE domain_dims_module,    only: get_domain_dims              ! get nx, ny, nt, and nbands
  USE domain_decomp_module,  only: get_domain_decomp_indices    ! get MPI domain decomposition indices
  
  USE time_windows_module,   only: get_time_windows             ! get info on the rolling time windows
  USE time_windows_module,   only: export_time_to_multiforce    ! populate legacy multiforce modules

  USE read_spatial_attrs,    only: read_spatial_attributes      ! read lat/lon and cell overlap
  USE read_elevbands_module, only: read_elevbands               ! read elevation bands

  USE read_hydromet_module,  only: read_hydromet_metadata       ! read hydromet metadata

  USE alloc_domain_module,   only: allocate_domain_data         ! allocate space for data arrays in the domain structure
  USE alloc_domain_module,   only: set_legacy_arrays            ! copy arrays in the domain%data structure to legacy arrays 

  USE init_mizuRoute,        only: init_mizuroute_domain        ! initialize mizuRoute structures used by FUSE

  implicit none
  
  ! input
  type(fuse_info)     , intent(inout)               :: info            ! domain info
  type(domain_data)   , intent(inout)               :: domain          ! domain data
  
  ! output
  integer(i4b)        , intent(out)                 :: ierr            ! error code
  character(len=1024) , intent(out)                 :: message         ! error message
  
  ! ----- internal -----------------------------------------------------------------------
  CHARACTER(LEN=1024)                               :: CMESSAGE        ! error message
  ! --------------------------------------------------------------------------------------
  ierr=0; message='setup_domain/'

  ! ----- set paths and file names -------------------------------------------------------
  
  ! read fuse control file (set paths/filenames etc.)
  call read_fuse_control_file(info, ierr, cmessage) 
  if (ierr/=0)then; message=trim(message)//trim(cmessage); ierr=20; return; endif

  ! ----- read domain metadata -----------------------------------------------------------
  
  ! populate domain structure with dimension lengths
  !   -- nx_global, ny_global, nt_global, n_bands
  call get_domain_dims(info, ierr, cmessage)
  if (ierr/=0)then; message=trim(message)//trim(cmessage); ierr=20; return; endif

  ! ----- read grid info and define indices for MPI domain decomposition ------------------

  ! open NetCDF hydromet file
  ierr = nf90_open(trim(info%files%fname_netcdf_hmet), nf90_nowrite, info%files%ncid_hydromet)
  if (ierr/=0)then; message=trim(message)//' nf90_open failed: '//trim(nf90_strerror(ierr)); return; endif
  if(isPrint) print *, 'Open hydromet file:', trim(info%files%fname_netcdf_hmet)
  if(isPrint) print *, 'NCID_HYDROMET is', info%files%ncid_hydromet
 
  ! ----- Compute time indices for sim/eval windows and subperiod chunk size --------------
  !
  ! Reads the hydromet-file NetCDF time coordinate (and units), and:
  !   - builds a Julian-day time axis and timestep in days
  !   - maps the user-specified simulation/evaluation date ranges into index windows
  !     (and optional subperiod chunks) stored in info%time
  call get_time_windows(info%files%ncid_hydromet, info, ierr, cmessage)
  if (ierr/=0)then; message=trim(message)//trim(cmessage); ierr=20; return; endif

  ! ----- initialize the mizuRoute data structures used by FUSE --------------------------

  call init_mizuroute_domain(info, domain, ierr, cmessage)
  if (ierr/=0)then; message=trim(message)//trim(cmessage); ierr=20; return; endif

  ! ----- MPI decomposition of the spatial domain ----------------------------------------

  ! get indices for MPI decomposition of the spatial domain: y_start_global, ny_local 
  ! NOTE: These indices will be used later to read different subsets of hydromet data for different ranks
  call get_domain_decomp_indices(info)

  ! ----- Allocate space for domain data --------------------------------------------------

  ! allocate space for the arrays in the domain%data structure
  call allocate_domain_data(info, domain, ierr, cmessage)
  if (ierr/=0)then; message=trim(message)//trim(cmessage); ierr=20; return; endif
 
  ! ----- Read static spatial information -------------------------------------------------

  ! read lat/lon and store in the domain%data%coords structure
  call read_spatial_attributes(info%files%ncid_hydromet, info, domain, ierr, message)
  if (ierr/=0)then; message=trim(message)//trim(cmessage); ierr=20; return; endif

  ! read elevation bands information and store in the domain%data structure
  call read_elevbands(info, domain, ierr, cmessage)
  if (ierr/=0)then; message=trim(message)//trim(cmessage); ierr=20; return; endif

  ! ----- Initialize time-varying hydromet inputs -----------------------------------------

  call read_hydromet_metadata(info%files%ncid_hydromet, info, ierr, cmessage)
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
