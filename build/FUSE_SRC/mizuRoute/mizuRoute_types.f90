module mizuroute_types

  ! shared general datatypes (FUSE and mizuRoute)
  use nrtype,   only: i4b, lgt

  ! mizuRoute data types
  use dataTypes, only: mizu_var_dlength => var_dlength
  use dataTypes, only: mizu_var_ilength => var_ilength
  use dataTypes, only: mizu_var_clength => var_clength

  implicit none
  private

  public :: mizuroute_topology

  type :: mizuroute_topology

    ! Network dimensions
    integer(i4b) :: n_hru = 0
    integer(i4b) :: n_seg = 0

    ! River-network data
    type(mizu_var_dlength), allocatable :: hru(:)
    type(mizu_var_dlength), allocatable :: seg(:)
    type(mizu_var_ilength), allocatable :: hru2seg(:)
    type(mizu_var_ilength), allocatable :: ntopo(:)
    type(mizu_var_clength), allocatable :: pfaf(:)

    logical(lgt) :: is_initialized = .false.

  end type mizuroute_topology

end module mizuroute_types
