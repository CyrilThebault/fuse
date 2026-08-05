module mizuroute_types

  use nrtype, only: i4b, lgt

  use dataTypes, only: mizu_var_dlength => var_dlength
  use dataTypes, only: mizu_var_ilength => var_ilength
  use dataTypes, only: mizu_var_clength => var_clength
  use dataTypes, only: mizu_remap       => remap
  use dataTypes, only: mizu_runoff      => runoff

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

    type(mizuroute_topology) :: topology  ! static network topology and attributes
    type(mizu_runoff)        :: runoff    ! FUSE runoff in mizuRouts structures

  end type river_network_data

  !---------------------------------------------------------------------
  ! Spatial mappings between model discretizations
  !---------------------------------------------------------------------
  type :: spatial_remap_data

    type(mizu_remap) :: forcing
    type(mizu_remap) :: routing

  end type spatial_remap_data

end module mizuroute_types
