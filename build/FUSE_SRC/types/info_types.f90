module info_types

 use nrtype

 use multiparam_types, only: par_id

 private
 public :: cli_options
 public :: time_info
 public :: space_info
 public :: fuse_info

 ! --------------------------------------------------------------------------------------
 
 type :: mpi_info
   logical(lgt) :: enabled = .false.
   integer(i4b) :: rank    = 0
   integer(i4b) :: nproc   = 1
 end type mpi_info

 ! -------------------------------------------------------------------------------------

 ! options for the command-line interface

 type :: cli_options
   character(len=:), allocatable :: tag              ! string to add to output file
   character(len=:), allocatable :: control_file
   character(len=:), allocatable :: domain_id
   character(len=:), allocatable :: runmode          ! def/idx/opt/sce
   character(len=:), allocatable :: sets_file        ! for idx,opt
   integer(i4b)                  :: indx = -1        ! for idx
   character(len=:), allocatable :: restart_freq     ! y/m/d/e/never
   logical(lgt)                  :: show_version = .false.
   logical(lgt)                  :: show_help    = .false.
   character(len=:), allocatable :: param_name(:)    ! list of parameter names
   real(sp), allocatable         :: param_value(:)   ! list of parameter values
 end type cli_options

 ! -------------------------------------------------------------------------------------

 type :: space_info
 
   ! global dimensions (full forcing file)
   integer(i4b) :: nx_global = 1
   integer(i4b) :: ny_global = 1

   ! local dimensions (after MPI split)
   integer(i4b) :: nx_local = 1
   integer(i4b) :: ny_local = 1

   ! decomposition along y dimension
   integer(i4b) :: y_start_global = 1
   integer(i4b) :: y_end_global   = 1

   ! mode flag
   logical(lgt) :: grid_flag = .false.
 
 end type space_info

 ! -------------------------------------------------------------------------------------
 
 type :: time_info
 
   ! forcing axis (global)
   integer(i4b) :: nt_global = 0
  
   ! simulation & evaluation indices into forcing time axis
   integer(i4b) :: sim_beg  = 1
   integer(i4b) :: sim_end  = 1
   integer(i4b) :: eval_beg = 1
   integer(i4b) :: eval_end = 1
  
   ! derived lengths
   integer(i4b) :: nt_sim   = 0
  
   ! subperiod / windowing
   logical(lgt) :: use_subperiods = .false.
   integer(i4b) :: nt_window      = 0       ! (= numtim_sub)
   integer(i4b) :: nt_window_cur  = 0       ! runtime: current window length
  
   ! bookkeeping for time axis
   character(len=:), allocatable :: units
   real(sp)       :: jdate_ref = 0._sp
   real(sp), allocatable :: time_steps(:)   ! time since reference time (transferred to output)
   real(sp), allocatable :: jdate(:)        ! julian day for each forcing record
 
 end type time_info

 ! -------------------------------------------------------------------------------------

 type :: snow_info
   integer(i4b) :: n_bands = 0
 end type snow_info

 ! -------------------------------------------------------------------------------------
 
 type :: file_info

   ! directories
   character(len=:), allocatable :: setngs_path
   character(len=:), allocatable :: input_path
   character(len=:), allocatable :: output_path
  
   ! settings filenames (relative or absolute)
   character(len=:), allocatable :: forcinginfo
   character(len=:), allocatable :: constraints
   character(len=:), allocatable :: mod_numerix
   character(len=:), allocatable :: m_decisions
  
   ! domain-derived input suffixes
   character(len=:), allocatable :: suffix_forcing
   character(len=:), allocatable :: suffix_elev_bands
  
   ! actual input filenames for this domain (derived once dom_id known)
   character(len=:), allocatable :: forcing_file    ! dom_id//suffix_forcing
   character(len=:), allocatable :: elevbands_file  ! dom_id//suffix_elev_bands
  
   ! output base name + concrete outputs
   character(len=512) :: fname_tempry
   character(len=512) :: fname_netcdf_forc
   character(len=512) :: fname_netcdf_runs
   character(len=512) :: fname_netcdf_para

   ! NetCDF forcing file info
   integer(i4b)       :: ncid_forc = -9999  ! NetCDF file ID for forcing data

 end type file_info

 ! -------------------------------------------------------------------------------------

 type :: run_config

  ! provenance
  character(len=512) :: file_manager_file = ""

  ! CLI options
  type(cli_options)  :: cli_opts

  ! model selection
  character(len=:), allocatable :: fmodel_id

  ! model information
  integer(i4b)       :: nState = -9999
  integer(i4b)       :: nParam = -9999
 
  ! number of input variables (3 = ppt, temp, pet; 4 = + obsq)
  integer(i4b)       :: nInput

  ! list of model parameters
  type(par_id), allocatable  :: listParam(:)

  ! run flags
  logical(lgt) :: q_only = .false.

  ! requested time windows (strings as read from filemanager)
  character(len=:), allocatable :: date_start_sim
  character(len=:), allocatable :: date_end_sim
  character(len=:), allocatable :: date_start_eval
  character(len=:), allocatable :: date_end_eval
  character(len=:), allocatable :: numtim_sub_str

  ! parsed / derived values (optional convenience)
  integer(i4b) :: numtim_sub = -9999      ! parsed from numtim_sub_str

  ! output dimension for number of parameter sets
  integer(i4b) :: nSets

  ! calibration metrics and metric transformations
  character(len=:), allocatable :: metric
  character(len=:), allocatable :: transfo

  ! SCE settings (store as numeric types)
  integer(i4b)      :: maxn  = -9999
  integer(i4b)      :: kstop = -9999
  real(sp)          :: pcento = -9999._sp

  ! store raw strings too
  character(len=20) :: maxn_str  = ""
  character(len=20) :: kstop_str = ""
  character(len=20) :: pcento_str = ""

 end type run_config

 ! -------------------------------------------------------------------------------------
 ! -------------------------------------------------------------------------------------
 
 type :: fuse_info
   type(mpi_info)   :: mpi
   type(space_info) :: space
   type(time_info)  :: time
   type(snow_info)  :: snow
   type(file_info)  :: files
   type(run_config) :: config
 end type fuse_info

end module info_types
