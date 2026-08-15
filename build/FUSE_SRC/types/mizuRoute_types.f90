module mizuroute_types

  use nrtype, only: wp, i4b, lgt

  use dataTypes, only: mizu_var_dlength => var_dlength
  use dataTypes, only: mizu_var_ilength => var_ilength
  use dataTypes, only: mizu_var_clength => var_clength
  use dataTypes, only: mizu_remap       => remap
  use dataTypes, only: mizu_runoff      => runoff

  use dataTypes, only: mizu_RCHPRP      => RCHPRP
  use dataTypes, only: mizu_RCHTOPO     => RCHTOPO

  use dataTypes, ONLY: mizu_STRFLX      => STRFLX
  use dataTypes, ONLY: mizu_STRSTA      => STRSTA

  implicit none
  private

  public :: mizuroute_topology
  public :: river_network_data
  public :: spatial_remap_data

  !---------------------------------------------------------------------
  ! River-network topology and attributes
  !---------------------------------------------------------------------
  type :: mizuroute_topology

    integer(i4b) :: n_hru = 0
    integer(i4b) :: n_seg = 0

    type(mizu_var_dlength), allocatable :: hru(:)
    type(mizu_var_dlength), allocatable :: seg(:)
    type(mizu_var_ilength), allocatable :: hru2seg(:)
    type(mizu_var_ilength), allocatable :: ntopo(:)
    type(mizu_var_clength), allocatable :: pfaf(:)

    logical(lgt) :: is_initialized = .false.

  end type mizuroute_topology

  !---------------------------------------------------------------------
  ! Persistent river-network data
  !---------------------------------------------------------------------
  type :: river_network_data

    type(mizuroute_topology)        :: topology  ! static network topology and attributes
    type(mizu_runoff)               :: runoff    ! FUSE runoff in mizuRoute structures

    ! mizuRoute routing: reach properties and network topology
    type(mizu_RCHPRP),  allocatable :: param(:)  ! reach properties
    type(mizu_RCHTOPO), allocatable :: ntopo(:)  ! network topology

    ! mizuRoute routing state and fluxes
    type(mizu_STRSTA),  allocatable :: state(:)  ! model states
    type(mizu_STRFLX),  allocatable :: flux(:)   ! model fluxes

    ! routing workspace
    real(wp),           allocatable :: reach_inflow(:)  ! lateral inflow to each reach [m3/s]

  end type river_network_data

  !---------------------------------------------------------------------
  ! Spatial mappings between model discretizations
  !---------------------------------------------------------------------
  type :: spatial_remap_data

    type(mizu_remap) :: forcing
    type(mizu_remap) :: routing

  end type spatial_remap_data

end module mizuroute_types
