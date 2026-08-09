module domain_dims_module
  use nrtype
  use info_types, only: fuse_info
  implicit none
  private
  public :: get_domain_dims

contains

  subroutine get_domain_dims(info, ierr, message)

  implicit none
  
  type(fuse_info), intent(inout) :: info
  integer(i4b),    intent(out)   :: ierr
  character(*),    intent(out)   :: message
  
  character(len=1024) :: hmet_file ! hydromet file
  character(len=1024) :: elev_file ! elev bands file

  character(len=1024) :: cmessage
  integer(i4b)        :: dimLen
  
  ierr = 0
  message = "get_domain_metadata/"

  ! get filenames
  hmet_file = trim(info%files%input_path)//trim(info%files%hydromet_file)
  elev_file = trim(info%files%input_path)//trim(info%files%elevbands_file)

  ! read forcing dimensions
  call read_forcing_dimensions(hmet_file, info, ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 
  ! read number of elevation bands
  call nc_get_dimlen_from_file(elev_file, "elevation_band", dimlen, ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  info%snow%n_bands = dimLen

  end subroutine get_domain_dims


  ! ----- utility routines --------------------------------------------------------------------------------------------

  ! ----- read forcing dimensions -------------------------------------------------------------------------------------

  subroutine read_forcing_dimensions(filepath, info, ierr, message)

  use nrtype
  use netcdf
  use info_types, only: fuse_info
  implicit none

  character(*),    intent(in)    :: filepath
  type(fuse_info), intent(inout) :: info
  integer(i4b),    intent(out)   :: ierr
  character(*),    intent(out)   :: message

  integer(i4b) :: ncid
  integer(i4b) :: varid
  integer(i4b) :: ndims
  integer(i4b) :: dimids(NF90_MAX_VAR_DIMS)
  character(len=NF90_MAX_NAME) :: dimname
  integer(i4b) :: idim,dimlen
  integer(i4b) :: time_varid
  
  associate(precip_name => info%files%precip_name, &
            grid_flag   => info%space%grid_flag,   &
            nx_global   => info%space%nx_global,   &
            ny_global   => info%space%ny_global,   &
            nt_global   => info%time%nt_global,    &
            nInput      => info%config%nInput)

  ierr=0; message="read_forcing_dimensions/"

  ! --- open NetCDF file for reading (nf90_nowrite) ---
  ierr = nf90_open(trim(filepath), nf90_nowrite, ncid)
  if(ierr /= nf90_noerr) then
    message = trim(message)//"nf90_open failed: "//trim(nf90_strerror(ierr))// &
              " [file="//trim(filepath)//"]"
    return
  endif

  ! --- get dimension lengths from precip variable shape ---
  ierr = nf90_inq_varid(ncid, trim(precip_name), varid)
  if(ierr /= nf90_noerr) then
    message = trim(message)//"cannot find var '"//trim(precip_name)//"': "//trim(nf90_strerror(ierr))
    return
  endif

  ierr = nf90_inquire_variable(ncid, varid, ndims=ndims, dimids=dimids)
  if(ierr /= nf90_noerr) then
    message = trim(message)//"inquire_variable failed: "//trim(nf90_strerror(ierr))
    return
  endif

  ! --- get the length of the time dimension (expect it is the last dimension) ---
  ierr = nf90_inquire_dimension(ncid, dimids(ndims), name=dimname, len=nt_global)
  if(ierr /= nf90_noerr) then
    message = trim(message)//"inquire_dimension failed: "//trim(nf90_strerror(ierr))
    return
  endif

  ! --- check that the last dimension is time ---
  if(trim(dimname) /= "time")then
    message=trim(message)//"FUSE expects (…, time) ordering; i.e., time dimension is last"
    ierr=20; return
  endif

  ! --- require rank 2 or 3 and require time last (already checked earlier) ---
  if (ndims /= 2 .and. ndims /= 3) then
    message = trim(message)//"expected forcing var rank 2 (spat,time) or 3 (spat,spat,time)"
    ierr = 20; return
  endif
 
  ! ndims == 2: enforce (hru,time) in  the feature order
  !   -> the only non-time dim is the "hru"/feature dimension
  if (ndims == 2) then

    ierr = nf90_inquire_dimension(ncid, dimids(1), len=ny_global)
    if(ierr /= nf90_noerr) then
      message = trim(message)//"inquire_dimension failed: "//trim(nf90_strerror(ierr))
      return
    endif
    nx_global = 1
  
  ! ndims == 3: enforce (x,y,time) in the file order
  !   -> can be (y,x,time) also since the spatial dimensions are general
  else

    ierr = nf90_inquire_dimension(ncid, dimids(1), len=nx_global)
    if(ierr /= nf90_noerr) then
      message = trim(message)//"inquire_dimension failed: "//trim(nf90_strerror(ierr))
      return
    endif
  
    ierr = nf90_inquire_dimension(ncid, dimids(2), len=ny_global)
    if(ierr /= nf90_noerr) then
      message = trim(message)//"inquire_dimension failed: "//trim(nf90_strerror(ierr))
      return
    endif

  endif  ! (ndims=3)

  ! define grid
  ! TODO: includes point list of catchments, but logic not implemented yet
  grid_flag = nx_global > 1 

  ! set the number of input variables (3 = ppt, temp, pet; 4 = + obsq)
  nInput = merge(3,4,grid_flag)

  ! --- close NetCDF file --- 
  ierr = nf90_close(ncid)
  if(ierr /= nf90_noerr) then
    message = trim(message)//"nf90_close failed: "//trim(nf90_strerror(ierr))// &
              " [file="//trim(filepath)//"]"
    return
  endif

  end associate
  end subroutine read_forcing_dimensions

  ! ----- get dimension length from file ------------------------------------------------------------------------------


  subroutine nc_get_dimlen_from_file(filepath, dimname, dimlen, ierr, message)

  use netcdf, only: nf90_open, nf90_close, nf90_nowrite, &
                    nf90_inq_dimid, nf90_inquire_dimension, &
                    nf90_strerror, nf90_noerr
  implicit none

  ! inputs
  character(*), intent(in)  :: filepath
  character(*), intent(in)  :: dimname

  ! outputs
  integer(i4b), intent(out) :: dimlen
  integer(i4b), intent(out) :: ierr
  character(*), intent(out) :: message

  ! locals
  integer(i4b) :: ncid, dimid

  ierr    = 0
  dimlen  = -1
  message = "nc_get_dimlen_from_file/"

  ! open NetCDF file for reading (nf90_nowrite)
  ierr = nf90_open(trim(filepath), nf90_nowrite, ncid)
  if(ierr /= nf90_noerr) then
    message = trim(message)//"nf90_open failed: "//trim(nf90_strerror(ierr))// &
              " [file="//trim(filepath)//"]"
    return
  endif

  ! get dimension ID
  ierr = nf90_inq_dimid(ncid, trim(dimname), dimid)
  if(ierr /= nf90_noerr) then
    message = trim(message)//"nf90_inq_dimid failed: "//trim(nf90_strerror(ierr))// &
              " [dim="//trim(dimname)//"]"
    return
  endif

  ! get dimension length
  ierr = nf90_inquire_dimension(ncid, dimid, len=dimlen)
  if(ierr /= nf90_noerr) then
    message = trim(message)//"nf90_inquire_dimension failed: "//trim(nf90_strerror(ierr))// &
              " [dim="//trim(dimname)//"]"
    return
  endif

  ! close
  ierr = nf90_close(ncid)
  if(ierr /= nf90_noerr) then
    message = trim(message)//"nf90_close failed: "//trim(nf90_strerror(ierr))
    return
  endif
  
  end subroutine nc_get_dimlen_from_file

end module domain_dims_module
