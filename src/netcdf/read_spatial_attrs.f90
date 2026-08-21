module read_spatial_attrs

use nrtype
use info_types,   only: fuse_info
use domain_types, only: coord_data
use domain_types, only: domain_data

use netcdf

use fuse_globaldata, only: NA_VALUE

implicit none

private

public :: read_spatial_attributes

contains

  subroutine read_spatial_attributes(ncid, info, domain, ierr, message)

  implicit none

  integer(i4b),      intent(in)    :: ncid
  type(fuse_info),   intent(in)    :: info
  type(domain_data), intent(inout) :: domain
  integer(i4b),      intent(out)   :: ierr
  character(*),      intent(out)   :: message

  character(len=256) :: cmessage

  ierr    = 0
  message = 'read_spatial_attributes/'

  ! latitude and longitude
  call read_latlon_2d(ncid, info, domain%coords, ierr, cmessage)
  if (ierr /= 0) then
    message = trim(message)//trim(cmessage)
    return
  endif

  ! area of each spatial element overlapping the basin
  call read_overlap_area(ncid, info, domain%olap_area, ierr, cmessage)
  if (ierr /= 0) then
    message = trim(message)//trim(cmessage)
    return
  endif

  end subroutine read_spatial_attributes

  ! --------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------

  subroutine read_latlon_2d(ncid, info, coord, ierr, message)

  implicit none

  integer(i4b),     intent(in)    :: ncid
  type(fuse_info),  intent(in)    :: info
  type(coord_data), intent(inout) :: coord
  integer(i4b),     intent(out)   :: ierr
  character(*),     intent(out)   :: message

  integer(i4b) :: vid_lat, vid_lon
  integer(i4b) :: nx, ny, ystart

  ierr    = 0
  message = 'read_latlon_2d/'

  nx     = info%space%nx_local
  ny     = info%space%ny_local
  ystart = info%space%y_start_global

  ! allocate local coordinate arrays
  allocate(coord%lat_2d(nx,ny), stat=ierr)
  if (ierr /= 0) then
    message = trim(message)//'unable to allocate latitude array'
    return
  endif

  allocate(coord%lon_2d(nx,ny), stat=ierr)
  if (ierr /= 0) then
    message = trim(message)//'unable to allocate longitude array'
    return
  endif

  ! get coordinate variable IDs
  ierr = nf90_inq_varid(ncid, trim(info%files%latitude_name), vid_lat)
  if (ierr /= nf90_noerr) then
    message = trim(message)//'unable to find latitude variable'
    return
  endif

  ierr = nf90_inq_varid(ncid, trim(info%files%longitude_name), vid_lon)
  if (ierr /= nf90_noerr) then
    message = trim(message)//'unable to find longitude variable'
    return
  endif

  ! read according to spatial layout already determined in get_domain_dims
  
  !                           is_gridded    is_clinear
  !  HRU / point list           F               F
  !  regular lat/lon grid       T               F
  !  rotated/curvilinear grid   T               T

  if (.not. info%space%is_gridded) then

    call read_latlon_point(ncid, vid_lat, vid_lon, &
                           ystart, ny,              &
                           coord%lat_2d, coord%lon_2d, ierr)

  elseif (info%space%is_clinear) then

    call read_latlon_curv(ncid, vid_lat, vid_lon, &
                          nx, ny, ystart,          &
                          coord%lat_2d, coord%lon_2d, ierr)

  else

    call read_latlon_rect(ncid, vid_lat, vid_lon, &
                          nx, ny, ystart,          &
                          coord%lat_2d, coord%lon_2d, ierr)

  endif

  if (ierr /= nf90_noerr) then
    message = trim(message)//'unable to read latitude/longitude: '// &
              trim(nf90_strerror(ierr))
    return
  endif

  end subroutine read_latlon_2d

  ! --------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------

  subroutine read_overlap_area(ncid, info, olap_area, ierr, message)

  implicit none

  integer(i4b),    intent(in)    :: ncid
  type(fuse_info), intent(in)    :: info
  real(wp),        intent(inout) :: olap_area(:,:)
  integer(i4b),    intent(out)   :: ierr
  character(*),    intent(out)   :: message

  integer(i4b) :: varid
  integer(i4b) :: nx, ny, ystart

  ierr    = 0
  message = 'read_overlap_area/'

  nx     = info%space%nx_local
  ny     = info%space%ny_local
  ystart = info%space%y_start_global

  ! get variable ID
  ierr = nf90_inq_varid(ncid, 'cell_area_in_basin', varid)
  if (ierr /= nf90_noerr) then
    message = trim(message)//"unable to find variable 'cell_area_in_basin'"
    return
  endif

  ! read local data
  if (info%space%is_gridded) then

    ierr = nf90_get_var(ncid, varid, olap_area, &
                        start=(/1,ystart/), count=(/nx,ny/))

  else

    ierr = nf90_get_var(ncid, varid, olap_area(1,:), &
                        start=(/ystart/), count=(/ny/))

  endif

  if (ierr /= nf90_noerr) then
    message = trim(message)//'unable to read cell_area_in_basin: '// &
              trim(nf90_strerror(ierr))
    return
  endif

  end subroutine read_overlap_area

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
  ! --------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------

  subroutine read_latlon_point(ncid, vid_lat, vid_lon, ystart, ny, &
                               lat_2d, lon_2d, ierr)

  implicit none

  integer(i4b), intent(in)  :: ncid, vid_lat, vid_lon
  integer(i4b), intent(in)  :: ystart, ny
  real(wp),     intent(out) :: lat_2d(:,:), lon_2d(:,:)
  integer(i4b), intent(out) :: ierr

  ierr = nf90_get_var(ncid, vid_lat, lat_2d(1,:), &
                      start=(/ystart/), count=(/ny/))
  if (ierr /= nf90_noerr) return

  ierr = nf90_get_var(ncid, vid_lon, lon_2d(1,:), &
                      start=(/ystart/), count=(/ny/))

  end subroutine read_latlon_point

  ! --------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------

  subroutine read_latlon_rect(ncid, vid_lat, vid_lon, nx, ny, ystart, &
                              lat_2d, lon_2d, ierr)

  implicit none

  integer(i4b), intent(in)  :: ncid, vid_lat, vid_lon
  integer(i4b), intent(in)  :: nx, ny, ystart
  real(wp),     intent(out) :: lat_2d(:,:), lon_2d(:,:)
  integer(i4b), intent(out) :: ierr

  real(wp) :: lat(ny)
  real(wp) :: lon(nx)

  ierr = nf90_get_var(ncid, vid_lat, lat, &
                      start=(/ystart/), count=(/ny/))
  if (ierr /= nf90_noerr) return

  ierr = nf90_get_var(ncid, vid_lon, lon, &
                      start=(/1/), count=(/nx/))
  if (ierr /= nf90_noerr) return

  lat_2d = spread(lat, dim=1, ncopies=nx)
  lon_2d = spread(lon, dim=2, ncopies=ny)

  end subroutine read_latlon_rect

  ! --------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------

  subroutine read_latlon_curv(ncid, vid_lat, vid_lon, nx, ny, ystart, &
                              lat_2d, lon_2d, ierr)
  
  implicit none

  integer(i4b), intent(in)  :: ncid, vid_lat, vid_lon
  integer(i4b), intent(in)  :: nx, ny, ystart
  real(wp),     intent(out) :: lat_2d(:,:), lon_2d(:,:)
  integer(i4b), intent(out) :: ierr

  ierr = nf90_get_var(ncid, vid_lat, lat_2d, &
                      start=(/1,ystart/), count=(/nx,ny/))
  if (ierr /= nf90_noerr) return

  ierr = nf90_get_var(ncid, vid_lon, lon_2d, &
                      start=(/1,ystart/), count=(/nx,ny/))

  end subroutine read_latlon_curv

  ! --------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------

end module read_spatial_attrs
