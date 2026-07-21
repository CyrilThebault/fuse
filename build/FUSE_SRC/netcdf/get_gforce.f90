module get_gforce_module

use nrtype
use info_types, only: fuse_info

use netcdf
use time_io

implicit none

private

public::get_varid
public::get_gforce_3d
public::read_latlon_2d

contains

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
  ierr = nf90_inq_varid(ncid, "latitude", vid_lat)
  if(ierr /= nf90_noerr) then
    message = trim(message)//"missing var 'latitude': "//trim(nf90_strerror(ierr))
    return
  endif

  ierr = nf90_inq_varid(ncid, "longitude", vid_lon)
  if(ierr /= nf90_noerr) then
    message = trim(message)//"missing var 'longitude': "//trim(nf90_strerror(ierr))
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
  
  SUBROUTINE get_varID(ncid,ierr,message)
  ! ---------------------------------------------------------------------------------------
  ! Creator:
  ! --------
  ! Nans Addor, 2017
  ! ---------------------------------------------------------------------------------------
  ! Purpose:
  ! --------
  ! Get NetCDF ID for each variable of the forcing file
  ! ---------------------------------------------------------------------------------------
  ! Modules Modified:
  ! -----------------
  ! MODULE multiforce -- populate structure ncid_var%(*)
  ! ---------------------------------------------------------------------------------------
  USE multiforce, only: nForce, nInput                   ! number of forcing variables
  USE multiforce, only: ncid_var                         ! NetCDF forcing variable ID
  
  USE multiforce,only:forcefile                          ! name of forcing file
  USE multiforce,only:vname_aprecip                      ! variable name: precipitation
  USE multiforce,only:vname_airtemp                      ! variable name: temperature
  USE multiforce,only:vname_spechum                      ! variable name: specific humidity
  USE multiforce,only:vname_airpres                      ! variable name: surface pressure
  USE multiforce,only:vname_swdown                       ! variable name: downward shortwave radiation
  USE multiforce,only:vname_potevap                      ! variable name: potential ET
  USE multiforce,only:vname_q                            ! variable indice: observed discharge
  
  USE multiforce,only:ilook_aprecip                      ! variable indice: precipitation
  USE multiforce,only:ilook_airtemp                      ! variable indice: temperature
  USE multiforce,only:ilook_spechum                      ! variable indice: specific humidity
  USE multiforce,only:ilook_airpres                      ! variable indice: surface pressure
  USE multiforce,only:ilook_swdown                       ! variable indice: downward shortwave radiation
  USE multiforce,only:ilook_potevap                      ! variable indice: potential ET
  USE multiforce,only:ilook_q                            ! variable indice: observed discharge
  
  IMPLICIT NONE
  
  ! input
  integer(i4b), intent(in)                :: ncid        ! NetCDF file ID
  ! output
  integer(i4b), intent(out)               :: ierr        ! error code
  character(*), intent(out)               :: message     ! error message
  ! internal
  integer(i4b),parameter                  :: strLen=1024 ! length of character string
  type names
   character(len=strLen)                  :: vname       ! singlecharacter strings
  end type names
  type(names),dimension(nForce)           :: cVec        ! names of character strings
  integer(i4b)                            :: iVar        ! loop through forcing data
  
  ! ---------------------------------------------------------------------------------------
  ! initialize error control
  ierr=0; message='get_varID/'
  
  ! get the vector of variable names
  cVec(ilook_aprecip)%vname = trim(vname_aprecip)  ! variable name: precipitation
  cVec(ilook_potevap)%vname = trim(vname_potevap)  ! variable name: potential ET
  cVec(ilook_airtemp)%vname = trim(vname_airtemp)  ! variable name: temperature
  cVec(ilook_q)%vname = trim(vname_q)              ! variable name: observed discharge
  cVec(ilook_spechum)%vname = trim(vname_spechum)  ! variable name: specific humidity
  cVec(ilook_airpres)%vname = trim(vname_airpres)  ! variable name: surface pressure
  cVec(ilook_swdown)%vname  = trim(vname_swdown)   ! variable name: downward shortwave radiation
  
  do ivar=1,nInput
 
    ! get the variable ID
    ierr = nf90_inq_varid(ncid, trim(cVec(iVar)%vname), ncid_var(ivar))
    if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'[variable='//trim(cVec(iVar)%vname)//']'; return; endif
  
  END DO
  
  END SUBROUTINE get_varID
  
  SUBROUTINE get_gforce_3d(info, itim_start, numtim, ierr, message)
  ! ---------------------------------------------------------------------------------------
  ! Creator:
  ! --------
  ! Nans Addor, based on Martyn Clark's get_gforce
  ! ---------------------------------------------------------------------------------------
  ! Purpose:
  ! --------
  ! Read NetCDF gridded forcing data for a range of time steps
  ! ---------------------------------------------------------------------------------------
  ! Modules Modified:
  ! -----------------
  ! MODULE multiforce -- populate structure GFORCE_3d(*,*)%(*)
  ! ---------------------------------------------------------------------------------------
  USE fuse_fileManager,only:INPUT_PATH                   ! defines data directory
  USE multiforce,only:forcefile                          ! name of forcing file
  USE multiforce,only:nSpat1, nSpat2, numtim_sub         ! dimensions of local slice
  USE multiforce,only:startSpat2                         ! starting y index for data read
  USE multiforce,only:vname_aprecip                      ! variable name: precipitation
  USE multiforce,only:vname_airtemp                      ! variable name: temperature
  USE multiforce,only:vname_spechum                      ! variable name: specific humidity
  USE multiforce,only:vname_airpres                      ! variable name: surface pressure
  USE multiforce,only:vname_swdown                       ! variable name: downward shortwave radiation
  USE multiforce,only:vname_potevap                      ! variable name: potential ET
  USE multiforce,only:vname_q                            ! variable name: observed discharge
  
  USE multiforce,only:ilook_aprecip                      ! variable indice: precipitation
  USE multiforce,only:ilook_airtemp                      ! variable indice: temperature
  USE multiforce,only:ilook_spechum                      ! variable indice: specific humidity
  USE multiforce,only:ilook_airpres                      ! variable indice: surface pressure
  USE multiforce,only:ilook_swdown                       ! variable indice: downward shortwave radiation
  USE multiforce,only:ilook_potevap                      ! variable indice: potential ET
  USE multiforce,only:ilook_q                            ! variable indice: observed discharge
  
  USE multiforce,only:ncid_var                           ! NetCDF ID for forcing variables
  USE multiforce,only:amult_ppt,amult_pet                ! multipliers o convert to mm/day
  USE multiforce,only:gForce_3d                          ! gridded forcing data
  USE multiforce,only:ancilF_3d                          ! ancillary forcing data
  USE multiforce,only:nForce, nInput                     ! number of forcing variables
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
  real(sp),parameter                      :: amiss=-9999._sp ! value for missing data
  integer(i4b),parameter                  :: strLen=1024 ! length of character string
  integer(i4b)                            :: iVar        ! loop through forcing data
  real(sp),dimension(:,:,:),allocatable   :: gTemp       ! temporary 3d grid
  type names
   character(len=strLen)                  :: vname       ! singlecharacter strings
  end type names
  type(names),dimension(nForce)           :: cVec        ! names of character strings
  logical(lgt),dimension(nForce)          :: lCheck      ! check the existence of variables
 
  !     integer(i4b) :: nx, ny, ystart
  !     integer(i4b) :: start_3d(3), count_3d(3)
  !    
  !     nx     = info%space%nx_local
  !     ny     = info%space%ny_local
  !     ystart = info%space%y_start_global
  !    
  !     start_3d = (/ 1, ystart, itim_start/)
  !     count_3d = (/nx,     ny,     numtim/)


  ! ---------------------------------------------------------------------------------------
  ! initialize error control
  ierr=0; message='get_gforce_3d/'
  ! ---------------------------------------------------------------------------------------
  
  ! initialize lCheck
  lCheck=.false.
 
  ! allocate space for the temporary grid
  allocate(gTemp(nSpat1,nSpat2,numtim), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating space for gTemp'; return; endif
  
  ! get the vector of variable names
  cVec(ilook_aprecip)%vname = trim(vname_aprecip)  ! variable name: precipitation
  cVec(ilook_potevap)%vname = trim(vname_potevap)  ! variable name: potential ET
  cVec(ilook_airtemp)%vname = trim(vname_airtemp)  ! variable name: temperature
  cVec(ilook_q)%vname = trim(vname_q)              ! variable name: observed discharge
  cVec(ilook_spechum)%vname = trim(vname_spechum)  ! variable name: specific humidity
  cVec(ilook_airpres)%vname = trim(vname_airpres)  ! variable name: surface pressure
  cVec(ilook_swdown)%vname  = trim(vname_swdown)   ! variable name: downward shortwave radiation
  
  ! get forcing grids
  do ivar=1,nInput
  
   ! get the data
   ierr = nf90_get_var(info%files%ncid_forc, ncid_var(ivar), gTemp, start=(/1,startSpat2,itim_start/), count=(/nSpat1,nSpat2,numtim/)); CALL HANDLE_ERR(IERR)
   if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

   ! save the data in the structure -- and convert fluxes to mm/day
   if(trim(cVec(iVar)%vname) == trim(vname_aprecip) )then
     
     gForce_3d(:,:,1:numtim)%ppt = gTemp(:,:,:)*amult_ppt; lCheck(ilook_aprecip) = .true.
  
   endif
  
   if(trim(cVec(iVar)%vname) == trim(vname_potevap) )then
       gForce_3d(:,:,1:numtim)%pet = gTemp(:,:,:)*amult_pet; lCheck(ilook_potevap) = .true.
   endif
  
   if(trim(cVec(iVar)%vname) == trim(vname_airtemp) )then
     gForce_3d(:,:,1:numtim)%temp = gTemp(:,:,:);       lCheck(ilook_airtemp) = .true.
   endif
  
   if(trim(cVec(iVar)%vname) == trim(vname_q) )then
     aValid(:,:,1:numtim)%obsq = gTemp(:,:,:);       lCheck(ilook_q) = .true.

   endif
  
   ! save the other variables required to compute PET
   !if( trim(cVec(iVar)%vname) == trim(vname_airtemp) )then; ancilF(:,:)%airtemp = gTemp(:,:,1);       lCheck(ilook_airtemp) = .true.; endif
   !if( trim(cVec(iVar)%vname) == trim(vname_spechum) )then; ancilF(:,:)%spechum = gTemp(:,:,1);       lCheck(ilook_spechum) = .true.; endif
   !if( trim(cVec(iVar)%vname) == trim(vname_airpres) )then; ancilF(:,:)%airpres = gTemp(:,:,1);       lCheck(ilook_airpres) = .true.; endif
   !if( trim(cVec(iVar)%vname) == trim(vname_swdown)  )then; ancilF(:,:)%swdown  = gTemp(:,:,1);       lCheck(ilook_swdown)  = .true.; endif
  
  end do  ! (loop thru forcing variables)
 
  ! deallocate space for gTemp
  deallocate(gTemp, stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem deallocating space for gTemp'; return; endif
  
   !PRINT *, 'PET', gForce_3d(:,:,1:numtim)%pet
   !PRINT *, 'PPT', gForce_3d(:,:,1:numtim)%ppt
   !PRINT *, 'TEMP', gForce_3d(:,:,1:numtim)%temp
  
  end subroutine get_gforce_3d

end module get_gforce_module
