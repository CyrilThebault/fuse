module network_routing_module

  use nrtype

  use dataTypes,       only: mizu_remap       => remap
  use dataTypes,       only: mizu_runoff      => runoff
  use mizuroute_types, only: mizuroute_topology

  use fuse_globaldata, only: isPrint
  use fuse_globaldata, only: do_mizuRoute

  implicit none

  private
  public :: network_routing

contains

  subroutine network_routing(runoff_data, routing_map, topology, ierr, message)

    use process_remap_module, only: remap_runoff

    implicit none

    ! input

    type(mizu_runoff)        , intent(inout)   :: runoff_data
    type(mizu_remap)         , intent(in)      :: routing_map
    type(mizuroute_topology) , intent(in)      :: topology

    ! output
    integer(i4b)             , intent(out)     :: ierr
    character(*)             , intent(out)     :: message

    ! locals
    integer(i4b) :: ihru
    integer(i4b) :: iseg
    character(len=256)  :: cmessage

    ! initialize error control
    ierr    = 0
    message = 'network_routing/'

    ! early return (not running mizuRoute)
    if ( .not. do_mizuRoute ) then
      if (isPrint) print*, 'mizuRoute hydrofabric file not defined: running lumped simulations'
      return
    endif

    ! remap runoff to basin HRUs
    call remap_runoff (runoff_data,             &   ! input: routed runoff from FUSE
                       routing_map,             &   ! input: mapping structure for routing
                       runoff_data%basinRunoff, &   ! output: runoff for basin HRUs
                       ierr, cmessage)              ! output: error control
    if (ierr /= 0) then; message = trim(message)//trim(cmessage); return; end if

  end subroutine network_routing

end module network_routing_module
