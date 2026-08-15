module network_routing_module

  use nrtype

  use dataTypes,       only: mizu_remap       => remap
  use dataTypes,       only: mizu_runoff      => runoff
  use mizuroute_types, only: mizuroute_topology
  use mizuroute_types, only: river_network_data
  use mizuroute_types, only: spatial_remap_data

  use fuse_globaldata, only: isPrint
  use fuse_globaldata, only: do_mizuRoute
  use fuse_globaldata, only: do_remapping

  implicit none

  private
  public :: network_routing

contains

  subroutine network_routing(river_network, routing_map, ierr, message)

    use process_remap_module, only: remap_runoff
    use process_remap_module, only: basin2reach

    implicit none

    ! input

    type(river_network_data) , intent(inout)   :: river_network
    type(mizu_remap)         , intent(in)      :: routing_map

    ! output
    integer(i4b)             , intent(out)     :: ierr
    character(*)             , intent(out)     :: message

    ! locals
    integer(i4b) :: iSeg
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
    if (do_remapping) then
      call remap_runoff(river_network%runoff,             &   ! input: routed runoff from FUSE
                        routing_map,                      &   ! input: mapping structure for routing
                        river_network%runoff%basinRunoff, &   ! output: runoff for basin HRUs
                        ierr, cmessage)                       ! output: error control
      if (ierr /= 0) then; message = trim(message)//trim(cmessage); return; end if
    end if

    ! map the basin runoff to the stream network
    call basin2reach(river_network%runoff%basinRunoff,    & ! input: basin runoff (m/s)
                     river_network%ntopo,                 & ! input: reach topology
                     river_network%param,                 & ! input: reach parameter
                     river_network%reach_inflow,          & ! output: reach inflow (m3/s)
                     ierr, cmessage)                        ! output: error control
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    print*, 'river_network%reach_inflow = ', river_network%reach_inflow
    stop

    ! transfer lateral inflow (routing workspace) to flux data structures
    do iSeg = 1,size(river_network%reach_inflow)
      river_network%flux(iSeg)%BASIN_QR(0) = river_network%flux(iSeg)%BASIN_QR(1)        ! streamflow from previous step
      river_network%flux(iSeg)%BASIN_QR(1) = river_network%reach_inflow(iSeg)            ! streamflow (m3/s)
    end do



    ! 

  end subroutine network_routing

end module network_routing_module
