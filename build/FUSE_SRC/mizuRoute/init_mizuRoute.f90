MODULE init_mizuRoute

! data types
USE nrtype,    ONLY: i4b,dp,lgt,strLen
USE dataTypes, ONLY: var_ilength     ! integer type:          var(:)%dat
USE dataTypes, ONLY: var_clength     ! integer type:          var(:)%dat
USE dataTypes, ONLY: var_dlength     ! double precision type: var(:)%dat, or dat

! metadata on data structures
USE globalData, ONLY: meta_HRU       ! HRU properties
USE globalData, ONLY: meta_HRU2SEG   ! HRU-to-segment mapping
USE globalData, ONLY: meta_SEG       ! stream segment properties
USE globalData, ONLY: meta_NTOPO     ! network topology
USE globalData, ONLY: meta_PFAF      ! pfafstetter code

! indices of named variables
USE var_lookup, ONLY: ixHRU      , nVarsHRU
USE var_lookup, ONLY: ixHRU2SEG  , nVarsHRU2SEG
USE var_lookup, ONLY: ixSEG      , nVarsSEG
USE var_lookup, ONLY: ixNTOPO    , nVarsNTOPO
USE var_lookup, ONLY: ixPFAF     , nVarsPFAF

! Shared data
USE public_var, ONLY: iulog
USE public_var, ONLY: charMissing
USE public_var, ONLY: integerMissing
USE public_var, ONLY: realMissing

! FUSE global variables
USE fuse_globaldata, only: isPrint
USE fuse_globaldata, only: do_mizuRoute

implicit none

private
public :: init_mizuroute_domain

CONTAINS

 !-----------------------------------------------------------------------
 ! Initialize the mizuRoute data structures used by FUSE.
 !
 ! This routine:
 !   (1) initializes the mizuRoute metadata;
 !   (2) configures the FUSE–mizuRoute interface;
 !   (4) constructs the river-network topology;
 !   (3) reads the spatial remapping information; and
 !   (5) allocates the routing input data structures.
 !-----------------------------------------------------------------------
 subroutine init_mizuroute_domain(info, domain, ierr, message)

  ! data types
  use info_types, only: fuse_info
  use domain_types, only: domain_data

  ! shared data
  use public_var, only: ancil_dir
  use public_var, only: idSegOut
  use public_var, only: ntopAugmentMode
  use globalData, only: onRoute

  ! external subroutines
  use popMetadat_module,   only: popMetadat           ! populate metadata
  use read_param_module,   only: read_param           ! read the routing parameters
  use read_remap,          only: get_remap_data       ! read remap data

  use nr_utils,            only: match_index

  implicit none

  type(fuse_info),   intent(in)    :: info
  type(domain_data), intent(inout) :: domain
  integer(i4b),      intent(out)   :: ierr
  character(*),      intent(out)   :: message

  integer(i4b)                     :: nSpace(1:2) = integerMissing
  character(len=strLen)            :: cmessage

  integer(i4b)                     :: iHRU
  integer(i4b), allocatable        :: basinID(:)

  ierr = 0
  message = 'init_mizuroute_domain/'

  ! set flag to run mizuRoute
  do_mizuRoute = allocated(info%ntopo%hfabric_file)

  ! early return (not running mizuRoute)
  if ( .not. do_mizuRoute ) then
    if (isPrint) print*, 'mizuRoute hydrofabric file not defined: running lumped simulations'
    return
  endif

  ! shortcuts to data structures
  associate(topology => domain%river_network%topology, &
            remap    => domain%remap%routing)

  !---------------------------------------------------------------------
  ! Read the mizuRoute namelist
  !---------------------------------------------------------------------

  call read_param(trim(info%mrout%namelist_path)//trim(info%mrout%namelist_file), ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  !---------------------------------------------------------------------
  ! Initialize mizuRoute metadata
  !---------------------------------------------------------------------
  
  ! Populate the default metadata structures
  call popMetadat(ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  !---------------------------------------------------------------------
  ! Configure the FUSE–mizuRoute interface
  !---------------------------------------------------------------------

  ! get the spatial dimensions
  nSpace(1) = info%space%ny_global ! latitude dimension
  nSpace(2) = info%space%nx_global ! longitude dimension

  ! Populate the shared mizuRoute control variables.
  call populate_mizu_modules(info)

  !---------------------------------------------------------------------
  ! Construct the river network topology
  !---------------------------------------------------------------------

  ! Write an augmented hydrofabric if an output filename is provided.
  ntopAugmentMode = allocated(info%ntopo%hfabric_newfile) 

  ! Enable all mizuRoute routing formulations during network initialization so that
  ! the complete set of routing-specific network data structures is available.
  onRoute(:) = .true.

  ! Read the hydrofabric and compute the derived network attributes.

  ! NOTE: init_ntopo is copied directly from mizuRoute without modification.
  !
  ! It is the only substantial mizuRoute routine duplicated in the FUSE compatibility layer; all other
  ! mizuRoute functionality is called from the original mizuRoute modules and subroutines.

  call init_ntopo(topology%n_hru,           &
                  topology%n_seg,           &
                  topology%hru,             &
                  topology%seg,             &
                  topology%hru2seg,         &
                  topology%ntopo,           &
                  topology%pfaf,            &
                  ierr, cmessage)

  if (ierr /= 0) then
    message = trim(message)//trim(cmessage)
    return
  end if

  topology%is_initialized = .true.

  !---------------------------------------------------------------------
  ! Read spatial remapping information
  !---------------------------------------------------------------------

  ! This defines the mapping between the FUSE hydrologic spatial units and the routing HRUs.
  
  if ( allocated(info%remap%remap_file) )then
   
    ! read runoff mapping file 
    call get_remap_data(trim(ancil_dir)//trim(info%remap%remap_file), & ! input: file name
                        nSpace,                                       & ! input: vector of spatial dimensions
                        remap,                                        & ! output: data structure to remap data from a polygon
                        ierr, cmessage)                                 ! output: error control
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! extract the vector of HRU IDs from the data structure
    basinID = [ (topology%hru2seg(iHRU)%var(ixHRU2SEG%hruId)%dat(1), iHRU=1,size(topology%hru2seg)) ]

    ! get indices of the HRU ids in the mapping file in the routing layer
    remap%hru_ix = match_index(basinID, remap%hru_id, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  endif  ! (if remapping file exists)

  end associate

  !---------------------------------------------------------------------
  ! Initialize routing input data structures
  !---------------------------------------------------------------------

  associate(runoff => domain%river_network%runoff, &
            n_hru  => domain%river_network%topology%n_hru)

  runoff%nSpace    = nSpace
  runoff%fillvalue = realMissing
  
  ! 1-D HRU runoff
  if ( .not. info%space%grid_flag ) then
    message=trim(message)//'HRU spatial config not yet implemented'
    ierr=10; return

  ! 2-D gridded runoff
  else
    allocate(runoff%sim2d(nSpace(1), nSpace(2)), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'unable to allocate gridded runoff input'; return; endif
  endif

  ! allocate space for HRU variables
  allocate(runoff%basinRunoff(n_hru), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'unable to allocate hru runoff input'; return; endif

  end associate

 end subroutine init_mizuroute_domain

 ! ---------------------------------------------------------------------
 ! ---------------------------------------------------------------------
 ! ---------------------------------------------------------------------
 ! ---------------------------------------------------------------------
 ! ---------------------------------------------------------------------
 ! ---------------------------------------------------------------------

 ! *********************************************************************
 ! private subroutine: provide information expected in mizuRoute modules
 ! *********************************************************************
 subroutine populate_mizu_modules(info)
 
  use info_types, only: fuse_info

  ! mizuRoute configuration expected by the unmodified source code
  !
  ! File paths/names
  use public_var, only: ancil_dir
  use public_var, only: fname_ntopOld
  use public_var, only: fname_ntopNew

  ! dimension names in hydrofabric file
  use public_var, only: dname_sseg
  use public_var, only: dname_nhru

  ! dimension names in remapping file
  use public_var, only: dname_hru_remap       ! name of dimension of river network HRU ID
  use public_var, only: dname_data_remap      ! name of dimension of runoff HRU overlapping with river network HRU
  
  ! variable names in remapping file
  use public_var, only: vname_hruid_in_remap  ! name of variable containing ID of river network HRU
  use public_var, only: vname_weight          ! name of variable contating areal weights of runoff HRUs within each river network HRU
  use public_var, only: vname_num_qhru        ! name of variable containing numbers of runoff HRUs within each river network HRU
  use public_var, only: vname_i_index         ! name of variable containing index of xlon dimension in runoff grid (if runoff file is grid)
  use public_var, only: vname_j_index         ! name of variable containing index of ylat dimension in runoff grid (if runoff file is grid)

  ! Routing options
  use public_var, only: idSegOut
  use public_var, only: ntopAugmentMode

  ! time step
  use public_var, only: dt

  implicit none

  type(fuse_info),   intent(in)    :: info

  ! -------------------------------------------------------------------
  ! Copy the hydrofabric settings read from the FUSE control file
  ! into the variables expected by the unmodified mizuRoute routines.
  ! -------------------------------------------------------------------

  ! path/name of hydrofabric file
  ancil_dir     = trim(info%ntopo%hfabric_path)
  fname_ntopOld = trim(info%ntopo%hfabric_file)

  ! name of augmented hydrofabric file 
  if(ntopAugmentMode) fname_ntopNew = trim(info%ntopo%hfabric_newfile)

  ! names of dimensions in the NetCDF file
  dname_sseg = trim(info%ntopo%dname_seg)
  dname_nhru = trim(info%ntopo%dname_hru)

  ! VARIABLE NAMES for data (overwrite default name in popMeta.f90)
  ! HRU structure
  meta_HRU    (ixHRU%area            )%varName = trim(info%ntopo%varname_area)      ! HRU area
  ! Mapping from HRUs to stream segments
  meta_HRU2SEG(ixHRU2SEG%HRUid       )%varName = trim(info%ntopo%varname_HRUid)     ! HRU id
  meta_HRU2SEG(ixHRU2SEG%hruSegId    )%varName = trim(info%ntopo%varname_hruSegId)  ! the stream segment id below each HRU
  ! network topology
  meta_NTOPO  (ixNTOPO%segId         )%varName = trim(info%ntopo%varname_segId)     ! unique id of each stream segment
  meta_NTOPO  (ixNTOPO%downSegId     )%varName = trim(info%ntopo%varname_downSegId) ! unique id of the next downstream segment
  ! reach properties
  meta_SEG    (ixSEG%length          )%varName = trim(info%ntopo%varname_length)    ! length of segment  (m)
  meta_SEG    (ixSEG%slope           )%varName = trim(info%ntopo%varname_slope)     ! slope of segment   (-)

  ! DIMENSION NAMES for remapping (overwrite default name in public_var.f90)
  dname_hru_remap      = trim(info%remap%dname_hru)         ! dimension name for river network HRU
  dname_data_remap     = trim(info%remap%dname_data)        ! dimension name for runoff HRU ID
  
  ! VARIABLE NAMES for remapping (overwrite default name in public_var.f90)
  vname_hruid_in_remap = trim(info%remap%vname_hruid)       ! variable name for river network hru id
  vname_weight         = trim(info%remap%vname_weight)      ! variable name for areal weights of runoff HRUs within each river network
  vname_num_qhru       = trim(info%remap%vname_num_qhru)    ! variable for numbers of runoff HRUs within each river network HRU
  vname_i_index        = trim(info%remap%vname_i_index)     ! variable for numbers of y (latitude) index if runoff file is grid
  vname_j_index        = trim(info%remap%vname_j_index)     ! variable for numbers of x (longitude) index if runoff file is grid

  ! network topology
  idSegOut = info%ntopo%idSegOut

  ! time step
  dt = info%mrout%dt

 end subroutine populate_mizu_modules

 ! *********************************************************************
 ! private subroutine: initialize river network data
 ! *********************************************************************
 SUBROUTINE init_ntopo(nHRU_out, nRch_out,                                           & ! output: number of HRU and Reaches
                       structHRU, structSEG, structHRU2SEG, structNTOPO, structPFAF, & ! output: data structure for river data
                       ierr, message)                                                  ! output: error controls
  ! Shared data
  USE public_var, ONLY: ancil_dir                ! name of the ancillary directory
  USE public_var, ONLY: fname_ntopOld            ! name of the old network topology file
  USE public_var, ONLY: fname_ntopNew            ! name of the new network topology file
  USE public_var, ONLY: dname_nhru               ! dimension name for HRUs
  USE public_var, ONLY: dname_sseg               ! dimension name for stream segments
  USE public_var, ONLY: maxPfafLen               ! maximum digit of pfafstetter code (default 32)
  ! options
  USE public_var, ONLY: ntopAugmentMode          ! River network augmentation mode
  USE public_var, ONLY: idSegOut                 ! River network subset mode (idSegOut > 0)
  ! global data
  USE globalData, ONLY: meta_PFAF                ! meta for pfafstetter code
  ! external subroutines
  USE read_streamSeg,       ONLY: getData                  ! get the ancillary data
  USE write_streamSeg,      ONLY: writeData                ! write the ancillary data
  USE process_ntopo,        ONLY: check_river_properties   ! check if river network data is physically valid
  USE ncio_utils,           ONLY: get_var_dims
  USE process_ntopo,        ONLY: augment_ntopo            ! compute all the additional network topology (only compute option = on)

  implicit none
  ! Argument variables
  integer(i4b)                  , intent(out) :: nHRU_out                 ! number of HRUs
  integer(i4b)                  , intent(out) :: nRch_out                 ! number of reaches
  type(var_dlength), allocatable, intent(out) :: structHRU(:)             ! HRU properties
  type(var_dlength), allocatable, intent(out) :: structSeg(:)             ! stream segment properties
  type(var_ilength), allocatable, intent(out) :: structHRU2SEG(:)         ! HRU-to-segment mapping
  type(var_ilength), allocatable, intent(out) :: structNTOPO(:)           ! network topology
  type(var_clength), allocatable, intent(out) :: structPFAF(:)            ! pfafstetter code
  integer(i4b)      , intent(out)             :: ierr                     ! error code
  character(*)      , intent(out)             :: message                  ! error message
  ! Local variables
  integer(i4b)                                :: tot_upstream             ! total number of all of the upstream stream segments for all stream segments
  integer(i4b)                                :: tot_upseg                ! total number of immediate upstream segments for all  stream segments
  integer(i4b)                                :: tot_hru                  ! total number of all the upstream hrus for all stream segments
  integer(i4b)                                :: tot_uh                   ! total number of unit hydrograph from all the stream segments
  integer(i4b),      allocatable              :: ixHRU_desired(:)         ! indices of desired hrus
  integer(i4b),      allocatable              :: ixSeg_desired(:)         ! indices of desired reaches
  integer(i4b)                                :: dummy(2)                 ! dummy variable to hold dimension length for 2D variables in netCDF
  integer(i4b)   , parameter                  :: maxUpstreamFile=90000000 ! 90 million: maximum number of upstream reaches to enable writing
  character(len=strLen)                       :: cmessage                 ! error message of downwind routine

  ierr=0; message='init_ntopo/'

  ! get the variable dimensions
  ! NOTE: need to update maxPfafLen to the exact character size for pfaf code in netCDF
  if (meta_PFAF(ixPFAF%code)%varFile) then
    call get_var_dims(trim(ancil_dir)//trim(fname_ntopOld), & ! input: file name
                      trim(meta_PFAF(ixPFAF%code)%varName), & ! input: pfaf code variable name in netcdf
                      ierr, cmessage,                       & ! output: error control
                      dlen=dummy)                             ! output optional: dimension length
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
    maxPfafLen = dummy(1)
  end if

  call getData(trim(ancil_dir)//trim(fname_ntopOld), & ! input: file name
               dname_nhru,    & ! input: dimension name of the HRUs
               dname_sseg,    & ! input: dimension name of the stream segments
               nHRU_out,      & ! output: number of HRUs
               nRch_out,      & ! output: number of stream segments
               structHRU,     & ! output: ancillary data for HRUs
               structSeg,     & ! output: ancillary data for stream segments
               structHRU2seg, & ! output: ancillary data for mapping hru2basin
               structNTOPO,   & ! output: ancillary data for network topology
               structPFAF,    & ! output: ancillary data for pfafstetter code
               ierr,cmessage)   ! output: error control
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  call check_river_properties(structNTOPO, structHRU, structSEG, ierr, cmessage) ! input: data structure for physical river network data
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  call augment_ntopo(nHRU_out,                         & ! number of HRUs
                     nRch_out,                         & ! number of stream segments
                     structHRU,                        & ! ancillary data for HRUs
                     structSeg,                        & ! ancillary data for stream segments
                     structHRU2seg,                    & ! ancillary data for mapping hru2basin
                     structNTOPO,                      & ! ancillary data for network toopology
                     ierr, cmessage,                   & ! error control
                     tot_hru       = tot_hru,          & ! total number of all the upstream hrus for all stream segments
                     tot_upseg     = tot_upseg,        & ! total number of all the immediate upstream segments for all stream segments
                     tot_upstream  = tot_upstream,     & ! total number of all the upstream segments for all stream segments
                     tot_uh        = tot_uh,           & ! total number of unit hydrograph for all stream segments
                     ixHRU_desired = ixHRU_desired,    & ! indices of desired hrus
                     ixSeg_desired = ixSeg_desired)      ! indices of desired reaches
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! write network topology (if augment mode or subset mode)
  if(ntopAugmentMode .or. idSegOut>0)then

    ! disable the dimension containing all upstream reaches
    ! NOTE: For the CONUS this is 1,872,516,819 reaches !!
    !        --> it will always be quicker to recompute than read+write
    !        --> users can modify the hard-coded parameter "maxUpstreamFile" if desired
    if(tot_upstream > maxUpstreamFile) tot_upstream=0

    call writeData(trim(ancil_dir)//trim(fname_ntopNew), & ! input: file name
                   tot_hru,       & ! input: total number of all the upstream hrus for all stream segments
                   tot_upseg,     & ! input: total number of immediate upstream segments for all  stream segments
                   tot_upstream,  & ! input: total number of all of the upstream stream segments for all stream segments
                   tot_uh,        & ! input: total number of unit hydrograph for all stream segments
                   ixHRU_desired, & ! input: indices of desired hrus
                   ixSeg_desired, & ! input: indices of desired reaches
                   structHRU,     & ! input: ancillary data for HRUs
                   structSeg,     & ! input: ancillary data for stream segments
                   structHRU2seg, & ! input: ancillary data for mapping hru2basin
                   structNTOPO,   & ! input: ancillary data for network topology
                   structPFAF,    & ! input: ancillary data for pfafstetter code
                   ierr,cmessage) ! output: error control
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    if (idSegOut>0) write(iulog,'(a)') 'Running in river network subset mode'
    if (ntopAugmentMode) write(iulog,'(a)') 'Running in river network augmentation mode'
    write(iulog,'(a)') 'Created a new network topology file '//trim(fname_ntopNew)
    write(iulog,'(a)') ' --> Run again using the new network topology file '
    return
  endif

 END SUBROUTINE init_ntopo

END MODULE init_mizuRoute
