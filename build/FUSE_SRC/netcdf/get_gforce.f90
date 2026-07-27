module get_gforce_module

use nrtype
use info_types, only: fuse_info

use netcdf

use fuse_globaldata, only: NVAR_FORC
use fuse_globaldata, only: iPRECIP, iTEMP, iPET, iQOBS

implicit none

private

public :: get_gforce_3d
public :: read_latlon_2d
public :: get_forcing_varids

contains

  ! --------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------
  
  subroutine read_latlon_2d(ncid, info, coord, ierr, message)

  use netcdf
  use nrtype
  use info_types, only: fuse_info
  use data_types, only: coord_data
  implicit none

  integer(i4b),    intent(in)    :: ncid
  type(fuse_info), intent(in)    :: info
  type(coord_data),intent(inout) :: coord
  integer(i4b),    intent(out)   :: ierr
  character(*),    intent(out)   :: message

  integer(i4b) :: vid_lat, vid_lon
  integer(i4b) :: nd_lat, nd_lon
  integer(i4b) :: dimids_lat(NF90_MAX_VAR_DIMS), dimids_lon(NF90_MAX_VAR_DIMS)
  integer(i4b) :: nx, ny, ystart
  integer(i4b) :: len_lat, len_lon
  integer(i4b) :: start2(2), count2(2)
  real(sp), allocatable :: lon_1d(:), lat_1d(:)

  ierr = 0
  message = "read_latlon_2d/"

  nx     = info%space%nx_local
  ny     = info%space%ny_local
  ystart = info%space%y_start_global

  ! Ensure 2D storage exists
  if (.not. allocated(coord%lat_2d)) allocate(coord%lat_2d(nx, ny))
  if (.not. allocated(coord%lon_2d)) allocate(coord%lon_2d(nx, ny))

  ! --- get varids ---
  ierr = nf90_inq_varid(ncid, trim(info%files%latitude_name), vid_lat)
  if(ierr /= nf90_noerr) then
    message = trim(message)//"missing var '"//trim(info%files%latitude_name)//"': "//trim(nf90_strerror(ierr))
    return
  endif

  ierr = nf90_inq_varid(ncid, trim(info%files%longitude_name), vid_lon)
  if(ierr /= nf90_noerr) then
    message = trim(message)//"missing var '"//trim(info%files%longitude_name)//"': "//trim(nf90_strerror(ierr))
    return
  endif

  ! --- ranks/dims ---
  ierr = nf90_inquire_variable(ncid, vid_lat, ndims=nd_lat, dimids=dimids_lat)
  if(ierr /= nf90_noerr) then
    message = trim(message)//"inquire latitude failed: "//trim(nf90_strerror(ierr))
    return
  endif

  ierr = nf90_inquire_variable(ncid, vid_lon, ndims=nd_lon, dimids=dimids_lon)
  if(ierr /= nf90_noerr) then
    message = trim(message)//"inquire longitude failed: "//trim(nf90_strerror(ierr))
    return
  endif

  !----------------------------------------------------------------------------
  ! Case A: Rectilinear OR point-list (lat 1D, lon 1D)
  !----------------------------------------------------------------------------
  if (nd_lat == 1 .and. nd_lon == 1) then

    ! Read full 1D vectors (easiest because slice depends on grid shape)
    ! NOTE: do MPI slice later
    ierr = nf90_inquire_dimension(ncid, dimids_lat(1), len=len_lat)
    if(ierr /= nf90_noerr) then
      message = trim(message)//"inquire lat dim failed: "//trim(nf90_strerror(ierr))
      return
    endif
    ierr = nf90_inquire_dimension(ncid, dimids_lon(1), len=len_lon)
    if(ierr /= nf90_noerr) then
      message = trim(message)//"inquire lon dim failed: "//trim(nf90_strerror(ierr))
      return
    endif

    allocate(lat_1d(len_lat))
    allocate(lon_1d(len_lon))

    ierr = nf90_get_var(ncid, vid_lat, lat_1d)
    if(ierr /= nf90_noerr) then
      message = trim(message)//"get_var(latitude) failed: "//trim(nf90_strerror(ierr))
      return
    endif

    ierr = nf90_get_var(ncid, vid_lon, lon_1d)
    if(ierr /= nf90_noerr) then
      message = trim(message)//"get_var(longitude) failed: "//trim(nf90_strerror(ierr))
      return
    endif

    coord%is_curvilinear = .false.
    coord%is_point_list  = (info%space%nx_global == 1)  ! our convention

    if (coord%is_point_list) then
      ! Point-list/HRU: lat(hru), lon(hru) -> store as (1,ny_local)
      if (nx /= 1) then
        message = trim(message)//"point-list detected but nx_local /= 1"
        ierr = 20; return
      endif
      coord%lat_2d(1,:) = lat_1d(ystart:ystart+ny-1)
      coord%lon_2d(1,:) = lon_1d(ystart:ystart+ny-1)

    else
      ! Rectilinear grid: lat(ny), lon(nx) -> broadcast to 2D
      
      ! lon_1d is global length nx_global
      coord%lon_2d(:,:) = spread(lon_1d(1:nx), dim=2, ncopies=ny)

      ! lat_1d is global length ny_global; take this rank's slice then replicate across x
      coord%lat_2d(:,:) = spread(lat_1d(ystart:ystart+ny-1), dim=1, ncopies=nx)
    
    endif

    deallocate(lat_1d, lon_1d)

    return
  endif

  !----------------------------------------------------------------------------
  ! Case B: Curvilinear (lat 2D, lon 2D)
  !----------------------------------------------------------------------------
  if (nd_lat == 2 .and. nd_lon == 2) then

    ! Read local slab in file order: (spat1,spat2) with y split along dim2

    start2 = (/ 1, ystart /)
    count2 = (/ nx, ny /)

    ierr = nf90_get_var(ncid, vid_lat, coord%lat_2d, start=start2, count=count2)
    if(ierr /= nf90_noerr) then
      message = trim(message)//"get_var(latitude 2D) failed: "//trim(nf90_strerror(ierr))
      return
    endif

    ierr = nf90_get_var(ncid, vid_lon, coord%lon_2d, start=start2, count=count2)
    if(ierr /= nf90_noerr) then
      message = trim(message)//"get_var(longitude 2D) failed: "//trim(nf90_strerror(ierr))
      return
    endif

    coord%is_curvilinear = .true.
    coord%is_point_list  = .false.

    return
  endif

  !----------------------------------------------------------------------------
  ! Anything else is unsupported under preprocessing + layout rules
  !----------------------------------------------------------------------------
  ierr = 20
  write(message,'(a,i0,a,i0,a)') trim(message)// &
    "unsupported lat/lon ranks (lat_ndims=", nd_lat, ", lon_ndims=", nd_lon, &
    "). If coords include time, preprocess to remove time from latitude/longitude."

  end subroutine read_latlon_2d

  ! --------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------
  subroutine get_dimIds(ncid, varid, nexpect, varDimIDs, ierr, message)
  ! used to get the vector of dimension ids for a given variable
  
  implicit none
  
  ! input
  integer(i4b),intent(in)   :: ncid     ! NetCDF file ID
  integer(i4b),intent(in)   :: varid    ! NetCDF variable ID
  integer(i4b),intent(in)   :: nexpect  ! number of dimensions expected
  
  ! output
  integer(i4b),intent(out)  :: varDimIDs(nexpect)  ! vector of dimension IDs
  integer(i4b),intent(out)  :: ierr     ! error code
  character(*), intent(out) :: message  ! error message
  
  ! internal variables
  integer(i4b)              :: nVarDims ! number of dimensions for given variable
  
  ! initialize error control
  ierr=0; message='get_dimIds/'
  
  ! get number of dimensions
  ierr = nf90_inquire_variable(ncid, varid, ndims=nVarDims)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
  
  ! check number of dimensions
  if(nVarDims/=nexpect)then; message=trim(message)//'unexpected number of dimensions for variable'; return; endif
  
  ! get vector of dimension IDs
  ierr = nf90_inquire_variable(ncid, varid, dimids=varDimIDs(:nVarDims))
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
  
  end subroutine get_dimIds
  
  ! --------------------------------------------------------------------------------------
 
  subroutine get_forcing_varids(ncid, info, ierr, message)
  implicit none
  integer(i4b), intent(in)       :: ncid
  type(fuse_info), intent(inout) :: info
  integer(i4b), intent(out)      :: ierr
  character(*), intent(out)      :: message

  integer(i4b) :: ivar

  ierr = 0
  message = "get_forcing_varids/"

  ! get table of name/varid pairs (names set in TOML read)
  info%files%forc%name(iPRECIP) = info%files%precip_name
  info%files%forc%name(iTEMP)   = info%files%temp_name
  info%files%forc%name(iPET)    = info%files%pet_name
  info%files%forc%name(iQOBS)   = info%files%qobs_name

  info%files%forc%varid(:) = -1

  ! get varid for each forcing variable
  do ivar = 1, NVAR_FORC

    if(info%space%grid_flag .and. ivar == iQOBS) cycle ! skips qobs if a grid

    call lookup_varid(ncid, trim(info%files%forc%name(ivar)), info%files%forc%varid(ivar), ierr, message)
    if(ierr/=0) return

  end do  ! ivar

  contains
    
    subroutine lookup_varid(ncid, vname, vid, ierr, message)
  
    integer(i4b), intent(in)    :: ncid
    character(len=*), intent(in) :: vname
    integer(i4b), intent(inout) :: vid
    integer(i4b), intent(inout) :: ierr
    character(*), intent(inout) :: message

    if (len_trim(vname) == 0) then
      ierr = 20
      message = trim(message)//"empty variable name"
      return
    end if

    ierr = nf90_inq_varid(ncid, trim(vname), vid)
    if (ierr /= 0) then
      message = trim(message)//trim(nf90_strerror(ierr))//" [var="//trim(vname)//"]"
    end if
    
   end subroutine lookup_varid
  
  end subroutine get_forcing_varids

  
  SUBROUTINE get_gforce_3d(info, itim_start, numtim, ierr, message)
  ! ---------------------------------------------------------------------------------------
  ! Creator:
  ! --------
  ! Nans Addor, based on Martyn Clark's get_gforce
  ! Modified by Martyn Clark to simplify and update to new data structures, 02/2026
  ! ---------------------------------------------------------------------------------------
  ! Purpose:
  ! --------
  ! Read NetCDF gridded forcing data for a range of time steps
  ! ---------------------------------------------------------------------------------------
  ! Modules Modified:
  ! -----------------
  ! MODULE multiforce -- populate structure GFORCE_3d(*,*)%(*)
  ! ---------------------------------------------------------------------------------------
  USE multiforce,only:gForce_3d                          ! gridded forcing data
  USE multiforce,only:aValid                             ! time series of lumped forcing/response data
  
  IMPLICIT NONE
  
  ! input
  type(fuse_info), intent(in)             :: info        ! info data structure that holds spatial indices
  integer(i4b),    intent(in)             :: itim_start  ! index of model time step - start of the period to extract
  integer(i4b),    intent(in)             :: numtim      ! number of model time steps to extract
  
  ! output
  integer(i4b),    intent(out)            :: ierr        ! error code
  character(*),    intent(out)            :: message     ! error message
  
  ! internal
  integer(i4b)                            :: iVar        ! loop through forcing data
  real(sp),dimension(:,:,:),allocatable   :: gTemp       ! temporary 3d grid
  integer(i4b)                            :: nx, ny       ! grid dimensions
  integer(i4b)                            :: ystart       ! start index iin input file (for MPI)
  integer(i4b)                            :: start_3d(3)  ! start indices in NetCDF file
  integer(i4b)                            :: count_3d(3)  ! count in NetCDF file
  
  ! initialize error control
  ierr=0; message='get_gforce_3d/'
  ! ---------------------------------------------------------------------------------------
  
  ! 3-d grid dimensions 
  nx     = info%space%nx_local
  ny     = info%space%ny_local
  ystart = info%space%y_start_global  ! start index in input file (for MPI)
  
  ! indices for NetCDF rea
  start_3d = (/ 1, ystart, itim_start/)
  count_3d = (/nx,     ny,     numtim/)

  ! allocate space for the temporary grid
  allocate(gTemp(nx,ny,numtim), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating space for gTemp'; return; endif
  
  ! get forcing grids
  do ivar = 1, NVAR_FORC
  
   if(info%space%grid_flag .and. ivar == iQOBS) cycle ! skips qobs if a grid
  
   ! get the data
   ierr = nf90_get_var(info%files%ncid_forc, info%files%forc%varid(ivar), gTemp, start=start_3d, count=count_3d)
   if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

   ! save the data in the structure -- and convert fluxes to mm/day
   select case(ivar)

    case (iPRECIP); gForce_3d(:,:,1:numtim)%ppt  = gTemp(:,:,:)
    case (iTEMP)  ; gForce_3d(:,:,1:numtim)%temp = gTemp(:,:,:)
    case (iPET)   ; gForce_3d(:,:,1:numtim)%pet  = gTemp(:,:,:)
    case (iQOBS)  ; aValid(   :,:,1:numtim)%obsq = gTemp(:,:,:)  ! TODO: check dimensions (works for nx=1, ny=1)
    case default
      message=trim(message)//'unable to identify forcing variable'
      ierr=10; return

   end select  ! identify forcing variable
  
  end do  ! (loop thru forcing variables)
 
  ! deallocate space for gTemp
  deallocate(gTemp, stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem deallocating space for gTemp'; return; endif
  
  end subroutine get_gforce_3d

end module get_gforce_module
