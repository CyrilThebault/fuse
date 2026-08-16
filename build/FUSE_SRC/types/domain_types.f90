module domain_types

 use nrtype

 use multiforce_types,  only: FDATA
 use multibands_types,  only: BANDS_INFO, BANDS_VAR
 use multistate_types,  only: STATEV
 use multi_flux_types,  only: FLUXES
 use multiroute_types,  only: RUNOFF

 use mizuroute_types,   only: river_network_data
 use mizuroute_types,   only: spatial_remap_data

 private
 public :: coord_data, domain_data

  ! ------------------------------------------------------------------------------------- 
  
  type :: coord_data

    logical(lgt) :: is_curvilinear = .false.   ! true if lat/lon are 2D
    logical(lgt) :: is_point_list  = .false.   ! true if nx=1 and lat/lon are 1D over ny
   
    ! 2D rectilinear OR point-list
    real(wp)    , allocatable :: lon_1d(:)   ! nx or ny depending on layout
    real(wp)    , allocatable :: lat_1d(:)
   
    ! 2D curvilinear
    real(wp)    , allocatable :: lon_2d(:,:) ! (nx_local, ny_local)
    real(wp)    , allocatable :: lat_2d(:,:)

    ! IDs for HRU and stream segments
    integer(i4b), allocatable :: hru_id(:)
    integer(i4b), allocatable :: seg_id(:)

    ! optional IDs (int is usually safest)
    integer(i4b), allocatable :: cell_id(:,:)  ! always stored locally as (nx_local, ny_local)

  end type coord_data

  ! -------------------------------------------------------------------------------------
  
  type :: domain_data

    ! coordinate information
    type(coord_data)             :: coords

    ! 3D forcing window
    type(FDATA), allocatable      :: force(:,:,:)       ! gForce_3d (nx_local, ny_local, nt_window)

    ! 3D state window
    type(STATEV), allocatable     :: state(:,:,:)       ! gState_3d (nx_local, ny_local, nt_window+1)

    ! 3D flux window
    type(FLUXES), allocatable     :: flux(:,:,:)        ! w_flux_3d (nx_local, ny_local, nt_window)

    ! 3D routing window
    type(RUNOFF), allocatable     :: route(:,:,:)       ! AROUTE_3d (nx_local, ny_local, nt_window)

    ! 2D elevation information
    logical(lgt), allocatable     :: elev_mask(:,:)     ! elev_mask (nx_local, ny_local)
    real(wp),     allocatable     :: z_forcing(:,:)     ! Z_FORCING_grid (nx_local, ny_local)

    ! 3D snow-band information
    type(BANDS_INFO), allocatable :: bands_info(:,:,:)  ! MBANDS_INFO_3d (nx_local, ny_local, n_bands)

    ! 4D snow-band state window
    type(BANDS_VAR), allocatable  :: bands_var(:,:,:,:) ! MBANDS_VAR_4d (nx_local, ny_local, n_bands, nt_window+1)

    ! Remapping and netework routing
    type(river_network_data)      :: river_network      ! river_network%topology, river_network%runoff, ...
    type(spatial_remap_data)      :: remap              ! remap%forcing, remap%routing

  end type domain_data

  ! -------------------------------------------------------------------------------------


end module domain_types
