MODULE init_mizuRoute

USE nrtype,    ONLY: i4b,dp,lgt,strLen

! fuse data types
USE info_types,   ONLY: fuse_info
USE domain_types, ONLY: domain_data

! mizuRoute data types
USE dataTypes, ONLY: var_ilength     ! integer type:          var(:)%dat
USE dataTypes, ONLY: var_clength     ! integer type:          var(:)%dat
USE dataTypes, ONLY: var_dlength     ! double precision type: var(:)%dat, or dat

! combined fuse+mizuRoute
USE mizuroute_types, ONLY: routing_time_data
USE mizuroute_types, ONLY: river_network_data
USE mizuroute_types, ONLY: spatial_remap_data

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

! Named variables for routing methods
USE public_var, ONLY: nRouteMethods           ! 6: number of routing methods available
USE public_var, ONLY: accumRunoff             ! 0: runoff accumulation over all the upstream reaches
USE public_var, ONLY: impulseResponseFunc     ! 1: impulse response function
USE public_var, ONLY: kinematicWaveTracking   ! 2: Lagrangian kinematic wave
USE public_var, ONLY: kinematicWave           ! 3: kinematic wave
USE public_var, ONLY: muskingumCunge          ! 4: muskingum-cunge
USE public_var, ONLY: diffusiveWave           ! 5: diffusiveWave

! indices in the active routing-method vector
USE globalData, ONLY: idxSUM,idxIRF,idxKWT, &
                       idxKW,idxMC, idxDW

! FUSE global variables
USE fuse_globaldata, only: isPrint
USE fuse_globaldata, only: do_remapping
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
 !   (5) allocates the mizuRoute routing data structures.
 !-----------------------------------------------------------------------
 subroutine init_mizuroute_domain(info, domain, ierr, message)

  ! shared data
  use public_var, only: ancil_dir
  use public_var, only: idSegOut
  use public_var, only: ntopAugmentMode
  use globalData, only: length_conv,time_conv
  
  use globalData, only: onRoute
  use globaldata, only: nRoutes
  use globaldata, only: routeMethods

  ! mizuRoute shim (unmodified mizuRoute code)
  use init_model_data_shim, only: init_ntopo
  use init_model_data_shim, only: init_route_method

  ! external mizuRoute subroutines
  use popMetadat_module,   only: popMetadat           ! populate metadata
  use read_param_module,   only: read_param           ! read the routing parameters
  use process_ntopo,       only: put_data_struct      ! copy data to the new structures 
  use read_remap,          only: get_remap_data       ! read remap data

  use nr_utils,            only: match_index

  implicit none

  type(fuse_info),   intent(inout) :: info
  type(domain_data), intent(inout) :: domain
  integer(i4b),      intent(out)   :: ierr
  character(*),      intent(out)   :: message

  integer(i4b)                     :: nSpace(1:2) = integerMissing
  character(len=strLen)            :: cmessage

  integer(i4b)                     :: iHRU
  integer(i4b)                     :: iSeg
  integer(i4b)                     :: idxRoute
  integer(i4b), allocatable        :: basinID(:)

  ierr = 0
  message = 'init_mizuroute_domain/'

  ! ---- early return (not running mizuRoute) ----
  if ( .not. do_mizuRoute ) then
    if (isPrint) print*, 'mizuRoute hydrofabric file not defined: running lumped simulations'
    return
  endif

  ! ---- initialize unit conversions (multipliers) ----
  length_conv = 1.0e-3_dp            ! FUSE runoff length: mm -> m
  time_conv   = 1._dp / 86400.0_dp   ! FUSE runoff time:   day-1 -> s-1

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
  call populate_mizu_modules(info, domain%river_network%time, ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! initialize polymorphic routing structures
  call init_route_method(ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

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

  call init_ntopo(domain%river_network%topology%n_hru,           &
                  domain%river_network%topology%n_seg,           &
                  domain%river_network%topology%hru,             &
                  domain%river_network%topology%seg,             &
                  domain%river_network%topology%hru2seg,         &
                  domain%river_network%topology%ntopo,           &
                  domain%river_network%topology%pfaf,            &
                  ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  
  ! TEMPORARY: Copy the data to the old mizuRoute data structures
  call put_data_struct(domain%river_network%topology%n_seg,      & 
                       domain%river_network%topology%seg,        &
                       domain%river_network%topology%ntopo,      &
                       domain%river_network%param,               &
                       domain%river_network%ntopo,               &
                       ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  info%space%n_hru = domain%river_network%topology%n_hru
  info%space%n_seg = domain%river_network%topology%n_seg

  domain%coords%hru_id = [ (domain%river_network%topology%hru2seg(iHRU)%var(ixHRU2SEG%hruId)%dat(1), iHRU=1,info%space%n_hru) ]
  domain%coords%seg_id = [ (domain%river_network%topology%ntopo  (iSeg)%var(ixNTOPO%segId  )%dat(1), iSeg=1,info%space%n_seg) ]

  domain%river_network%topology%is_initialized = .true.

  !---------------------------------------------------------------------
  ! Read spatial remapping information
  !---------------------------------------------------------------------

  ! This defines the mapping between the FUSE hydrologic spatial units and the routing HRUs.
  
  if ( do_remapping ) then
   
    ! read runoff mapping file 
    call get_remap_data(trim(ancil_dir)//trim(info%remap%remap_file), & ! input: file name
                        nSpace,                                       & ! input: vector of spatial dimensions
                        domain%remap%routing,                         & ! output: data structure to remap data from a polygon
                        ierr, cmessage)                                 ! output: error control
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
    
      domain%remap%routing%hru_ix = match_index(domain%coords%hru_id, domain%remap%routing%hru_id, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  endif  ! (if remapping file exists)

  !---------------------------------------------------------------------
  ! Allocate mizuRoute routing structures
  !---------------------------------------------------------------------

  ! only enable mizuRoute routing formulations for those desired
  ! NOTE: do this after reading network topology because want ntopo for all methods
  onRoute(:) = .false.
  onRoute(routeMethods) = .true.
   
  call allocate_mizuroute_domain(info,                                 &
                                 domain%river_network,                 &
                                 nSpace,                               &
                                 ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

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
 subroutine populate_mizu_modules(info, time, ierr, message)
 
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
  use globaldata, only: routeMethods
  use globaldata, only: nRoutes
  use globalData, only: nMolecule

  ! time step for routing model
  use public_var, only: secprday
  use public_var, only: dt_route => dt  ! seconds

  use nr_utils,   only: char2int        ! convert a character string to an integer vector

  implicit none

  type(fuse_info),          intent(in)     :: info
  type(routing_time_data),  intent(out)    :: time

  integer(i4b),             intent(out)    :: ierr
  character(*),             intent(out)    :: message

  real(dp)                                 :: dt_fuse   ! FUSE time step (seconds)
  integer(i4b)                             :: iRoute

  ierr = 0
  message = "populate_mizu_modules/"

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

  ! routing methods
  call char2int(trim(info%mrout%methods), routeMethods, invalid_value=0)
  nRoutes = size(routeMethods)

  do iRoute = 1, nRoutes
    select case(routeMethods(iRoute))
      case (accumRunoff);           idxSUM = iRoute
      case (kinematicWaveTracking); idxKWT = iRoute
      case (impulseResponseFunc);   idxIRF = iRoute
      case (muskingumCunge);        idxMC  = iRoute
      case (kinematicWave);         idxKW  = iRoute
      case (diffusiveWave);         idxDW  = iRoute
      case default
        message=trim(message)//'routOpt may include invalid digits; expect digits 1-5 in routOpt'; ierr=81; return
    end select
  end do

  ! indices in the active routing-method vector
  do iRoute = 1, nRoutes
    if (routeMethods(iRoute)==kinematicWave) then
      nMolecule%KW_ROUTE = 20
    else if (routeMethods(iRoute)==muskingumCunge) then
      nMolecule%MC_ROUTE = 2
    else if (routeMethods(iRoute)==diffusiveWave) then
      nMolecule%DW_ROUTE = 20
    end if
  end do

  ! the number of computational "molecules" in each reach
  do iRoute = 1, nRoutes
    if (routeMethods(iRoute)==kinematicWave) then
      nMolecule%KW_ROUTE = 20
    else if (routeMethods(iRoute)==muskingumCunge) then
      nMolecule%MC_ROUTE = 2
    else if (routeMethods(iRoute)==diffusiveWave) then
      nMolecule%DW_ROUTE = 20
    end if
  end do

  ! network topology
  idSegOut = info%ntopo%idSegOut

  ! time step
  dt_route = info%mrout%dt

  ! -------------------------------------------------------------------
  ! set up time step lengths
  ! -------------------------------------------------------------------

  dt_fuse = info%time%deltim_days * secprday

  if (dt_route > dt_fuse) then
    dt_route = dt_fuse
    print*, 'WARNING: dt_route > dt_fuse; setting dt_route = dt_fuse'
  end if

  time%n_sub  = ceiling(dt_fuse / dt_route)
  time%dt_sub = dt_fuse / real(time%n_sub, dp)

  dt_route = time%dt_sub

 end subroutine populate_mizu_modules

 ! *********************************************************************
 ! private subroutine: allocate space for the mizuRoute structures
 ! *********************************************************************
 subroutine allocate_mizuroute_domain(info, river_network, nSpace, &
                                      ierr, message)
 
   use globalData, only: onRoute
   use globalData, only: nMolecule
  
   use globaldata, only: routeMethods
   use globaldata, only: nRoutes

   implicit none
   
   type(fuse_info),              intent(in)    :: info
   type(river_network_data),     intent(inout) :: river_network
   integer(i4b),                 intent(in)    :: nSpace(2)
   integer(i4b),                 intent(out)   :: ierr
   character(*),                 intent(out)   :: message
   
   character(len=strLen)                       :: cmessage
   
   integer(i4b)                                :: iHRU, n_hru
   integer(i4b)                                :: iSeg, n_seg
   integer(i4b)                                :: n_time
   integer(i4b)                                :: idxRoute

   n_hru  = info%space%n_hru
   n_seg  = info%space%n_seg
   n_time = info%time%nt_window

   ! ---- allocate space for runoff inputs ----
   
   river_network%runoff%nSpace    = nSpace
   river_network%runoff%fillvalue = realMissing
   
   ! 1-D HRU runoff
   if ( .not. info%space%grid_flag ) then
     message=trim(message)//'HRU spatial config not yet implemented'
     ierr=10; return
   
   ! 2-D gridded runoff
   else
     allocate(river_network%runoff%sim2d(nSpace(1), nSpace(2)), stat=ierr)
     if(ierr/=0)then; message=trim(message)//'unable to allocate gridded runoff input'; return; endif
   endif
   
   ! allocate space for HRU variables
   allocate(river_network%runoff%basinRunoff(n_hru), stat=ierr)
   if(ierr/=0)then; message=trim(message)//'unable to allocate hru runoff input'; return; endif
   
   ! ---- initialize network states and fluxes ----
   
   ! allocate space for all segments in the river network
   allocate(river_network%flux(n_seg),  &
            river_network%state(n_seg), & 
            river_network%reach_inflow(n_seg), stat=ierr)
   if(ierr/=0)then; message=trim(message)//'unable to allocate river_network flux/state'; return; endif
   
   ! * loop through stream segments
   do iSeg = 1, n_seg
   
     ! allocate fluxes for the routing method vector in each stream segment
     allocate(river_network%flux(iSeg)%ROUTE(nRoutes), stat=ierr)
     if (ierr /= 0) then
       write(message,'(A,I0)') trim(message)//'unable to allocate river_network%flux%ROUTE for iSeg=', iSeg
       return
     end if
   
     ! * loop through ACTIVE routing methods
     do idxRoute = 1, nRoutes
     
       ! allocate states each routing method individually
   
       select case(routeMethods(idxRoute))
   
         case (kinematicWave)
           allocate(river_network%state(iSeg)%KW_ROUTE%molecule%Q(nMolecule%KW_ROUTE), &
                    source=0._dp, stat=ierr)
       
         case (muskingumCunge)
           allocate(river_network%state(iSeg)%MC_ROUTE%molecule%Q(nMolecule%MC_ROUTE), &
                    source=0._dp, stat=ierr)
   
         case (diffusiveWave)
           allocate(river_network%state(iSeg)%DW_ROUTE%molecule%Q(nMolecule%DW_ROUTE), &
                    source=0._dp, stat=ierr)
         
         case (accumRunoff, impulseResponseFunc, kinematicWaveTracking)
           write(message,'(A,I0,A,I0)') trim(message)//'routing method =', routeMethods(idxRoute), &
                                        'not implemented in FUSE: use standalone mizuRoute'
           ierr=10; return

         case default
           message=trim(message)//'unable to identify routing method'
           ierr=10; return
   
       end select
   
       if (ierr /= 0) then
         write(message,'(A,I0,A,I0)') trim(message)//'unable to allocate routing state for iSeg=', &
                                      iSeg, ', method =', routeMethods(idxRoute)
         return
       endif
   
       ! initialize common routing inputs
       river_network%flux(iSeg)%BASIN_QR(:)   = 0._dp
       river_network%flux(iSeg)%REACH_WM_FLUX = 0._dp
       river_network%flux(iSeg)%REACH_WM_VOL  = 0._dp

       ! method-specific routing fluxes
       river_network%flux(iSeg)%ROUTE(idxRoute)%REACH_VOL   = 0._dp
       river_network%flux(iSeg)%ROUTE(idxRoute)%REACH_Q     = 0._dp
       river_network%flux(iSeg)%ROUTE(idxRoute)%Qerror      = 0._dp
     
     end do  ! * loop through routing methods
  end do  ! * loop through stream segments
 
  ! ---- allocate space for routing outputs ---- 

  allocate(river_network%method(nRoutes), stat=ierr)
  if (ierr /= 0) then
    message = trim(message)//'unable to allocate routing method data'
    return
  end if
  
  ! * loop through ACTIVE routing methods
  do idxRoute=1,nRoutes
     
    allocate(river_network%method(idxRoute)%streamflow(n_seg, n_time), &
             source=0._dp, stat=ierr)
  
    if (ierr /= 0) then
      write(message,'(A,I0)') &
        trim(message)//'unable to allocate streamflow for routing method=', routeMethods(idxRoute)
      return
    end if
  
  end do

 end subroutine allocate_mizuroute_domain 

END MODULE init_mizuRoute
