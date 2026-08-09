! Edited by Brian Henn to include snow model, 7/2013
! Edited by Nans Addor to set simulation and evaluation periods, 11/2017
! Modified by Martyn Clark to populate domain structure, 12/2025
MODULE fuse_filemanager

  use nrtype
  use info_types, only: cli_options
  use info_types, only: fuse_info

  use fuse_globaldata, only: NVAR_HYDROMET
  use fuse_globaldata, only: iPRECIP, iTEMP, iPET, iQOBS

  implicit none
  private

  public :: read_fuse_control_file

  ! ----- all of these legacy globals are stored in the "info" data structure ---------------------

  ! expose legacy globals
  public :: SETNGS_PATH, INPUT_PATH, OUTPUT_PATH
  public :: suffix_hydromet, suffix_elev_bands
  public :: M_DECISIONS, CONSTRAINTS, MOD_NUMERIX, MBANDS_NC
  public :: FMODEL_ID, Q_ONLY_STR, Q_ONLY
  public :: date_start_sim, date_end_sim, date_start_eval, date_end_eval, numtim_sub_str
  public :: METRIC, TRANSFO
  public :: KSTOP_str, MAXN_str, PCENTO_str

  ! ------------------------------------------------------------------------------------------------

  ! FUSE-wide pathlength
  integer(i4b),parameter::fusePathLen=512

  ! defines the path for data files
  CHARACTER(LEN=fusePathLen)  :: SETNGS_PATH
  CHARACTER(LEN=fusePathLen)  :: INPUT_PATH
  CHARACTER(LEN=fusePathLen)  :: OUTPUT_PATH
  
  ! content of input directory
  CHARACTER(LEN=fusePathLen)  :: suffix_hydromet   ! suffix for hydromet file
  CHARACTER(LEN=fusePathLen)  :: suffix_elev_bands ! suffix for elevation band file
  
  ! content of settings directory
  CHARACTER(LEN=fusePathLen)  :: M_DECISIONS       ! definition of model decisions
  CHARACTER(LEN=fusePathLen)  :: CONSTRAINTS       ! definition of parameter constraints
  CHARACTER(LEN=fusePathLen)  :: MOD_NUMERIX       ! definition of numerical solution technique
  CHARACTER(LEN=fusePathLen)  :: MBANDS_NC         ! netcdf file defining the elevation bands
  
  ! content of output directory
  CHARACTER(LEN=64)           :: FMODEL_ID         ! string defining FUSE model
  CHARACTER(LEN=64)           :: Q_ONLY_STR        ! TRUE = restrict attention to simulated runoff
  LOGICAL                     :: Q_ONLY            ! .TRUE. = restrict attention to simulated runoff
  
  ! define simulation and evaluation periods
  CHARACTER(len=20)           :: date_start_sim    ! date start simulation
  CHARACTER(len=20)           :: date_end_sim      ! date end simulation
  CHARACTER(len=20)           :: date_start_eval   ! date start evaluation period
  CHARACTER(len=20)           :: date_end_eval     ! date end evaluation period
  CHARACTER(len=20)           :: numtim_sub_str    ! number of time steps of subperiod (will be kept in memory)

  ! evaluation metrics and transformation
  CHARACTER(len=20)           :: METRIC            ! metric chosen as objective function
  CHARACTER(len=20)           :: TRANSFO           ! streamflow transformation

  ! SCE parameters
  CHARACTER(len=20)           :: KSTOP_str         ! number of shuffling loops the value must change by PCENTO
  CHARACTER(len=20)           :: MAXN_str          ! maximum number of trials before optimization is terminated
  CHARACTER(len=20)           :: PCENTO_str        ! the percentage

contains

  ! -------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------

  subroutine read_fuse_control_file(fuseFileManagerIn, opts, info, err, message)
  use tomlf_all, only: toml_table, toml_error, toml_key, toml_value ! data types
  use tomlf_all, only: toml_load, get_value                         ! procedures
  
  ! Purpose: Reads FUSE control file (TOML version) AND populates info structure
  !
  implicit none

  ! dummies
  character(*),      intent(in)     :: fuseFileManagerIn
  type(cli_options), intent(in)     :: opts
  type(fuse_info),   intent(inout)  :: info
  integer(i4b),      intent(out)    :: err
  character(*),      intent(out)    :: message

  ! TOML table
  type(toml_table), allocatable :: tbl         ! root TOML table
  type(toml_table), pointer     :: subtable    ! sub-table for a given section
  type(toml_key),   allocatable :: sections(:) ! top-level sections
  type(toml_key),   allocatable :: keys(:)     ! sub-table keys
  type(toml_error), allocatable :: error

  ! locals
  integer(i4b) :: istat
  integer(i4b) :: n, i, j
  character(len=256) :: lookup

  ! create file paths
  character(len=256) :: dom_id, tag, run_mode

  err = 0
  message = "read_fuse_control_file/"

  ! ----- load the root TOML table -----
  call toml_load(tbl, trim(fuseFileManagerIn), error=error)
  if (allocated(error)) then
    message = "problem loading toml file['"//trim(fuseFileManagerIn)//"']: "//trim(error%message)
    err=10; return
  endif

  ! ----- get the top-level sections -----
  call tbl%get_keys(sections)
  if(.not.allocated(sections)) then
    message = "problem loading toml sections['"//trim(fuseFileManagerIn)//"']"
    err=10; return
  endif

  ! ----- loop through sections -----
  do i = 1, size(sections)

    ! ----- load the TOML sub-table for the current section -----
    call get_value(tbl, trim(sections(i)%key), subtable, requested=.false.) 
    if(.not.associated(subtable)) then
      message = "problem loading toml sub-sections['"//trim(fuseFileManagerIn)//"']:"//trim(sections(i)%key)
      err=10; return
    endif

    ! ----- get keys for a given section (sub-table) -----
    call subtable%get_keys(keys)
    
    ! ----- loop through the sub-table -----
    do j = 1, size(keys)


      ! ---------------------------------------------------------------------------------------------------------------
      ! ---------------------------------------------------------------------------------------------------------------

      ! ----- control file assignment block ---------------------------------------------------------------------------
      
      lookup = trim(sections(i)%key)//'.'//trim(keys(j)%key)

      select case(lookup)
        
        ! ---- files: paths ----
        case ("filepaths.input_dir"          ); call get_value(subtable, trim(keys(j)%key), info%files%input_path       , stat=istat)
        case ("filepaths.output_dir"         ); call get_value(subtable, trim(keys(j)%key), info%files%output_path      , stat=istat)
        case ("filepaths.settings_dir"       ); call get_value(subtable, trim(keys(j)%key), info%files%setngs_path      , stat=istat)

        ! ---- files: suffixes ----
        case ("input.hydromet_suffix"        ); call get_value(subtable, trim(keys(j)%key), info%files%suffix_hydromet  , stat=istat)
        case ("input.elevbands_suffix"       ); call get_value(subtable, trim(keys(j)%key), info%files%suffix_elev_bands, stat=istat)

        ! ---- files: settings filenames ----
        case ("model.decisions_file"         ); call get_value(subtable, trim(keys(j)%key), info%files%m_decisions      , stat=istat)
        case ("model.numerics_file"          ); call get_value(subtable, trim(keys(j)%key), info%files%mod_numerix      , stat=istat)
        case ("model.constraints_file"       ); call get_value(subtable, trim(keys(j)%key), info%files%constraints      , stat=istat)

        ! ---- files: forcing coordinate names ----
        case ("forcing_coords.time"          ); call get_value(subtable, trim(keys(j)%key), info%files%time_name        , stat=istat)
        case ("forcing_coords.latitude"      ); call get_value(subtable, trim(keys(j)%key), info%files%latitude_name    , stat=istat)
        case ("forcing_coords.longitude"     ); call get_value(subtable, trim(keys(j)%key), info%files%longitude_name   , stat=istat)

        ! ---- files: hydromet variable names ----
        case ("hydromet_vars.precip"         ); call get_value(subtable, trim(keys(j)%key), info%files%precip_name      , stat=istat)
        case ("hydromet_vars.temp"           ); call get_value(subtable, trim(keys(j)%key), info%files%temp_name        , stat=istat)
        case ("hydromet_vars.pet"            ); call get_value(subtable, trim(keys(j)%key), info%files%pet_name         , stat=istat)
        case ("hydromet_vars.qobs"           ); call get_value(subtable, trim(keys(j)%key), info%files%qobs_name        , stat=istat)

        ! ---- mizuRoute: namelist path/filenames ----
        case ("mizuRoute.namelist_path"      ); call get_value(subtable, trim(keys(j)%key), info%mrout%namelist_path    , stat=istat)
        case ("mizuRoute.namelist_file"      ); call get_value(subtable, trim(keys(j)%key), info%mrout%namelist_file    , stat=istat)

        ! ---- mizuRoute: runtime ----
        case ("mizuRoute.dt"                 ); call get_value(subtable, trim(keys(j)%key), info%mrout%dt               , stat=istat)
        
        ! ---- hydrofabric: path/filenames ----
        case ("hydrofabric.hfabric_path"     ); call get_value(subtable, trim(keys(j)%key), info%ntopo%hfabric_path     , stat=istat)
        case ("hydrofabric.hfabric_file"     ); call get_value(subtable, trim(keys(j)%key), info%ntopo%hfabric_file     , stat=istat)
        case ("hydrofabric.hfabric_newfile"  ); call get_value(subtable, trim(keys(j)%key), info%ntopo%hfabric_newfile  , stat=istat)

        ! ---- hydrofabric: dimensions ----
        case ("hydrofabric.dname_hru"        ); call get_value(subtable, trim(keys(j)%key), info%ntopo%dname_hru        , stat=istat)
        case ("hydrofabric.dname_seg"        ); call get_value(subtable, trim(keys(j)%key), info%ntopo%dname_seg        , stat=istat)

        ! ---- hydrofabric: variable names ----
        case ("hydrofabric.varname_HRUid"    ); call get_value(subtable, trim(keys(j)%key), info%ntopo%varname_HRUid    , stat=istat)
        case ("hydrofabric.varname_segId"    ); call get_value(subtable, trim(keys(j)%key), info%ntopo%varname_segId    , stat=istat)
        case ("hydrofabric.varname_hruSegId" ); call get_value(subtable, trim(keys(j)%key), info%ntopo%varname_hruSegId , stat=istat)
        case ("hydrofabric.varname_downSegId"); call get_value(subtable, trim(keys(j)%key), info%ntopo%varname_downSegId, stat=istat)
        case ("hydrofabric.varname_area"     ); call get_value(subtable, trim(keys(j)%key), info%ntopo%varname_area     , stat=istat)
        case ("hydrofabric.varname_slope"    ); call get_value(subtable, trim(keys(j)%key), info%ntopo%varname_slope    , stat=istat)
        case ("hydrofabric.varname_length"   ); call get_value(subtable, trim(keys(j)%key), info%ntopo%varname_length   , stat=istat)

        ! ---- hydrofabric: network topology ----
        case ("hydrofabric.seg_outlet"       ); call get_value(subtable, trim(keys(j)%key), info%ntopo%idSegOut         , stat=istat)

        ! ---- remapping: filename ----
        case ("remapping.remap_file"         ); call get_value(subtable, trim(keys(j)%key), info%remap%remap_file       , stat=istat)
 
        ! ---- remapping: dimension names ----
        case ("remapping.dname_hru"          ); call get_value(subtable, trim(keys(j)%key), info%remap%dname_hru        , stat=istat)
        case ("remapping.dname_data"         ); call get_value(subtable, trim(keys(j)%key), info%remap%dname_data       , stat=istat)
        
        ! ---- remapping: variable names ----
        case ("remapping.vname_hruid"        ); call get_value(subtable, trim(keys(j)%key), info%remap%vname_hruid      , stat=istat)
        case ("remapping.vname_weight"       ); call get_value(subtable, trim(keys(j)%key), info%remap%vname_weight     , stat=istat)
        case ("remapping.vname_num_qhru"     ); call get_value(subtable, trim(keys(j)%key), info%remap%vname_num_qhru   , stat=istat)
        case ("remapping.vname_i_index"      ); call get_value(subtable, trim(keys(j)%key), info%remap%vname_i_index    , stat=istat)
        case ("remapping.vname_j_index"      ); call get_value(subtable, trim(keys(j)%key), info%remap%vname_j_index    , stat=istat)

        ! ---- config: runtime ----
        case ("output.model_id"              ); call get_value(subtable, trim(keys(j)%key), info%config%fmodel_id       , stat=istat)
        case ("output.q_only"                ); call get_value(subtable, trim(keys(j)%key), info%config%q_only          , stat=istat)

        ! ---- config: periods ----
        case ("run_periods.date_start_sim"   ); call get_value(subtable, trim(keys(j)%key), info%config%date_start_sim  , stat=istat)
        case ("run_periods.date_end_sim"     ); call get_value(subtable, trim(keys(j)%key), info%config%date_end_sim    , stat=istat)
        case ("run_periods.date_start_eval"  ); call get_value(subtable, trim(keys(j)%key), info%config%date_start_eval , stat=istat)
        case ("run_periods.date_end_eval"    ); call get_value(subtable, trim(keys(j)%key), info%config%date_end_eval   , stat=istat)
        case ("run_periods.numtim_sub_str"   ); call get_value(subtable, trim(keys(j)%key), info%config%numtim_sub_str  , stat=istat)

        ! ---- config: calibration ----
        case ("calibration.metric"           ); call get_value(subtable, trim(keys(j)%key), info%config%metric          , stat=istat)
        case ("calibration.transfo"          ); call get_value(subtable, trim(keys(j)%key), info%config%transfo         , stat=istat)

        ! ---- config: SCE (read numeric, then next populate legacy strings) ----
        case ("sce.maxn"                     ); call get_value(subtable, trim(keys(j)%key), info%config%maxn            , stat=istat)
        case ("sce.kstop"                    ); call get_value(subtable, trim(keys(j)%key), info%config%kstop           , stat=istat)
        case ("sce.pcento"                   ); call get_value(subtable, trim(keys(j)%key), info%config%pcento          , stat=istat)

        ! ---- default case (something in the table that is not specified above) -----
        case default
          message = trim(message)// "unexpected entry: section = "//trim(sections(i)%key)//"; sub-section = "//trim(keys(j)%key)
          err=20; return
      
      end select ! (select key/value pair based on lookup)

      ! ---- error checking -----
      if(istat /= 0)then
        message=trim(message)// "get_value error: section = "//trim(sections(i)%key)//"; sub-section = "//trim(keys(j)%key)
        err=20; return
      endif

    end do  ! (looping through sub-sections)
  end do  ! (looping through sections)

  ! ---- populate legacy strings ----
  write(info%config%maxn_str,  '(i0)'     ) info%config%maxn
  write(info%config%kstop_str, '(i0)'     ) info%config%kstop
  write(info%config%pcento_str,'(es20.10)') info%config%pcento

  ! ---- domain id, run mode and tag for output files ----
  dom_id   = trim(opts%domain_id)
  run_mode = trim(opts%runmode)

  tag = ""
  if(allocated(opts%tag)) tag = trim(opts%tag)

  ! ---- derived input filenames ----
  info%files%hydromet_file  = trim(dom_id)//trim(info%files%suffix_hydromet)
  info%files%elevbands_file = trim(dom_id)//trim(info%files%suffix_elev_bands)

  ! ---- derived output base name ----
  info%files%fname_tempry = trim(info%files%output_path)// &
                            trim(dom_id)//'_'//trim(info%config%fmodel_id)//'_'//trim(tag)

  ! ---- final filenames ----
  info%files%fname_netcdf_hmet = trim(info%files%input_path)//trim(info%files%hydromet_file)
  info%files%fname_netcdf_runs = trim(info%files%fname_tempry)//'_runs_'//trim(run_mode)//'.nc'
  info%files%fname_netcdf_para = trim(info%files%fname_tempry)//'_para_'//trim(run_mode)//'.nc'

  ! ---- populate legacy modules ----
  call export_domain_to_legacy(info)

  end subroutine read_fuse_control_file

  ! -------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------
  
  ! ----- export domain config variables to legacy modules ------------------------------

  subroutine export_domain_to_legacy(info)
  use model_defn, only: FNAME_TEMPRY, FNAME_NETCDF_RUNS, FNAME_NETCDF_PARA
  use multiparam, only: MAXN, KSTOP, PCENTO

  implicit none

  type(fuse_info), intent(in)  :: info

  ! ---- populate legacy module globals  ----

  SETNGS_PATH       = trim(info%files%setngs_path)
  INPUT_PATH        = trim(info%files%input_path)
  OUTPUT_PATH       = trim(info%files%output_path)
  suffix_hydromet   = trim(info%files%suffix_hydromet)
  suffix_elev_bands = trim(info%files%suffix_elev_bands)
  CONSTRAINTS       = trim(info%files%constraints)
  MOD_NUMERIX       = trim(info%files%mod_numerix)
  M_DECISIONS       = trim(info%files%m_decisions)

  FMODEL_ID         = trim(info%config%fmodel_id)
  Q_ONLY            = info%config%q_only

  date_start_sim    = trim(info%config%date_start_sim)
  date_end_sim      = trim(info%config%date_end_sim)
  date_start_eval   = trim(info%config%date_start_eval)
  date_end_eval     = trim(info%config%date_end_eval)
  numtim_sub_str    = trim(info%config%numtim_sub_str)

  METRIC            = trim(info%config%metric)
  TRANSFO           = trim(info%config%transfo)

  MAXN              = info%config%maxn
  KSTOP             = info%config%kstop
  PCENTO            = info%config%pcento

  MAXN_str          = trim(info%config%maxn_str)
  KSTOP_str         = trim(info%config%kstop_str)
  PCENTO_str        = trim(info%config%pcento_str)

  ! populate module model_defn
  FNAME_TEMPRY      = trim(info%files%fname_tempry)
  FNAME_NETCDF_RUNS = trim(info%files%fname_netcdf_runs)
  FNAME_NETCDF_PARA = trim(info%files%fname_netcdf_para)
  
  ! populate shared public variable in this module (fuse_filemanager)
  MBANDS_NC = trim(info%files%elevbands_file)
  
  ! populate module multiparam
  MAXN   = info%config%maxn
  KSTOP  = info%config%kstop
  PCENTO = info%config%pcento
  
  ! ---- logging ----
  print *, 'Paths defined in file manager:'
  print *, 'SETNGS_PATH:', trim(info%files%setngs_path)
  print *, 'INPUT_PATH:',  trim(info%files%input_path)
  print *, 'OUTPUT_PATH:', trim(info%files%output_path)
  
  print *, 'Dates defined in file manager:'
  print *, 'date_start_sim:',  trim(info%config%date_start_sim)
  print *, 'date_end_sim:',    trim(info%config%date_end_sim)
  print *, 'date_start_eval:', trim(info%config%date_start_eval)
  print *, 'date_end_eval:',   trim(info%config%date_end_eval)
  print *, 'numtim_sub_str:',  trim(info%config%numtim_sub_str)
  
  print *, 'Q_ONLY', info%config%q_only
  
  end subroutine export_domain_to_legacy

!----------------------------------------------------
END MODULE fuse_filemanager
