MODULE init_mizuRoute_topo

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

! modules

implicit none

private
public :: init_mizuroute_topology

CONTAINS

 subroutine init_mizuroute_topology(info, domain, ierr, message)

  use info_types, only: fuse_info
  use data_types, only: domain_data
  
  use public_var, only: idSegOut
  use public_var, only: ntopAugmentMode
  use public_var, only: impulseResponseFunc
  use globalData, only: onRoute

  use  popMetadat_module,   only: popMetadat

  implicit none

  type(fuse_info),   intent(in)    :: info
  type(domain_data), intent(inout) :: domain
  integer(i4b),      intent(out)   :: ierr
  character(*),      intent(out)   :: message

  character(len=strLen) :: cmessage

  ierr = 0
  message = 'init_mizuroute_topology/'

  ! Populate the default metadata structures
  call popMetadat(ierr, cmessage)
  if (ierr /= 0) then
    message = trim(message)//trim(cmessage)
    return
  end if

  ! Populate the mizuRoute data modules
  call populate_mizu_modules(info)

  ! flag to write the augmented hydrofabric
  !  -- write if augmented filename provided in the control file
  ntopAugmentMode = allocated(info%ntopo%hfabric_newfile) 

  ! enable allocation of the network structures required for impulse-response-function routing.
  onRoute = .false.
  onRoute(impulseResponseFunc) = .true.

  ! Read the hydrofabric and compute the derived network attributes.
  ! NOTE: init_ntopo is copied direct from mizuRoute
  call init_ntopo(domain%river_network%n_hru,           &
                  domain%river_network%n_seg,           &
                  domain%river_network%hru,             &
                  domain%river_network%seg,             &
                  domain%river_network%hru2seg,         &
                  domain%river_network%ntopo,           &
                  domain%river_network%pfaf,            &
                  ierr, cmessage)

  if (ierr /= 0) then
    message = trim(message)//trim(cmessage)
    return
  end if

  domain%river_network%is_initialized = .true.

 end subroutine init_mizuroute_topology

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
  ! File paths/names and dimension/variable names
  use public_var, only: ancil_dir
  use public_var, only: fname_ntopOld
  use public_var, only: fname_ntopNew
  use public_var, only: dname_sseg
  use public_var, only: dname_nhru

  ! Routing options
  use public_var, only: idSegOut
  use public_var, only: ntopAugmentMode

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

  ! network topology
  idSegOut = info%ntopo%idSegOut

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

END MODULE init_mizuRoute_topo
