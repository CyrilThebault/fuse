module data_types

 use nrtype

 use multiforce_types,  only: ADATA, FDATA, VDATA
 use multibands_types,  only: BANDS_VAR
 use multistate_types,  only: STATEV
 use multi_flux_types,  only: FLUXES
 use multiroute_types,  only: RUNOFF

 private
 public :: coord_data, domain_data

  ! ------------------------------------------------------------------------------------- 
  
  type :: coord_data

    logical(lgt) :: is_curvilinear = .false.   ! true if lat/lon are 2D
    logical(lgt) :: is_point_list  = .false.   ! true if nx=1 and lat/lon are 1D over ny
   
    ! 2D rectilinear OR point-list
    real(sp), allocatable :: lon_1d(:)   ! nx or ny depending on layout
    real(sp), allocatable :: lat_1d(:)
   
    ! 2D curvilinear
    real(sp), allocatable :: lon_2d(:,:) ! (nx_local, ny_local)
    real(sp), allocatable :: lat_2d(:,:)
   
    ! optional IDs (int is usually safest)
    integer(i4b), allocatable :: cell_id(:,:)  ! always stored locally as (nx_local, ny_local)

  end type coord_data

  ! -------------------------------------------------------------------------------------
  
  type :: domain_data

    ! coordinate information
    type(coord_data)             :: coords

    ! 2D ancillary forcing (optional, for PET etc.)
    type(ADATA), allocatable     :: ancil(:,:)         ! (nx_local, ny_local)

    ! 3D forcing window (nx_local, ny_local, numtim_sub)
    type(FDATA), allocatable     :: force(:,:,:)       ! force_3d

    ! 3D state window (nx_local, ny_local, numtim_sub+1)
    type(STATEV), allocatable    :: state(:,:,:)       ! state_3d

    ! 3D flux window (nx_local, ny_local, numtim_sub)
    type(FLUXES), allocatable    :: flux(:,:,:)        ! flux_3d

    ! 3D routing window (nx_local, ny_local, numtim_sub)
    type(RUNOFF), allocatable    :: route(:,:,:)       ! route_3d

    ! 4D snow-band state window (nx_local, ny_local, n_bands, numtim_sub+1)
    type(BANDS_VAR), allocatable :: bands(:,:,:,:)      ! bands_var_4d

    ! 3D observed discharge / validity (optional)
    type(VDATA), allocatable     :: valid(:,:,:)       ! (nx_local, ny_local, numtim_sub)

    ! basin-average time series for output convenience
    type(FDATA), allocatable     :: aForce(:)          ! (numtim_sub)
    type(RUNOFF), allocatable    :: aRoute(:)          ! (numtim_sub)

  end type domain_data

  ! -------------------------------------------------------------------------------------
  
end module data_types
