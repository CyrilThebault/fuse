module read_elevbands_module
  !!
  !! Read elevation-band information (area fraction + mean elevation) from an elevation-band NetCDF file,
  !! supporting either:
  !!   - rank-3 variables: (spat1, spat2, band)   [grid]
  !!   - rank-2 variables: (spat2, band)          [list/HRU], mapped internally to (nx=1, ny=spat2, band)
  !!
  !! This module reads into REAL arrays (af, zmid) so you can then copy into your derived type
  !! MBANDS_INFO_3d(:,:,:)%AF / %Z_MID without baking that type into the reader.
  !!
  use nrtype
  use netcdf
  use info_types, only: fuse_info
  use info_types, only: space_info
  use data_types, only: domain_data
  implicit none
  private

  public :: read_elevbands

contains

  subroutine read_elevbands(info, domain, ierr, message)

    use fuse_globaldata,   only: NA_VALUE_SP

    implicit none

    type(fuse_info)     , intent(in)      :: info      ! domain info
    type(domain_data)   , intent(inout)   :: domain    ! domain data

    integer(i4b)        , intent(out)     :: ierr
    character(*)        , intent(out)     :: message

    character(len=1024) :: eb_file ! elev bands file

    integer(i4b) :: ncid_eb
    integer(i4b) :: vid_af, vid_me
    integer(i4b) :: nd_af, nd_me
    integer(i4b) :: i, j, ib
    real(wp) :: afsum
    real(wp) :: af(info%space%nx_local, info%space%ny_local, info%snow%n_bands)
    real(wp) :: zmid(info%space%nx_local, info%space%ny_local, info%snow%n_bands)
    character(len=1024) :: cmessage

    ierr = 0
    message = "read_elevbands_arrays/"

    ! ---- get name of file that holds elevation bands ----
    eb_file = trim(info%files%input_path)//trim(info%files%elevbands_file)

    ! --- open NetCDF file for reading (nf90_nowrite) ---
    ierr = nf90_open(trim(eb_file), nf90_nowrite, ncid_eb)
    if(ierr /= nf90_noerr) then
      message = trim(message)//"nf90_open failed: "//trim(nf90_strerror(ierr))// &
                " [file="//trim(eb_file)//"]"
      return
    endif

    ! ---- lookup varids ----
    ierr = nf90_inq_varid(ncid_eb, "area_frac", vid_af)
    if (ierr /= nf90_noerr) then
      message = trim(message)//"nf90_inq_varid(area_frac) failed: "//trim(nf90_strerror(ierr))
      return
    end if

    ierr = nf90_inq_varid(ncid_eb, "mean_elev", vid_me)
    if (ierr /= nf90_noerr) then
      message = trim(message)//"nf90_inq_varid(mean_elev) failed: "//trim(nf90_strerror(ierr))
      return
    end if

    ! ---- inquire ranks (we support 2D or 3D) ----
    ierr = nf90_inquire_variable(ncid_eb, vid_af, ndims=nd_af)
    if (ierr /= nf90_noerr) then
      message = trim(message)//"inquire_variable(area_frac) failed: "//trim(nf90_strerror(ierr))
      return
    end if

    ierr = nf90_inquire_variable(ncid_eb, vid_me, ndims=nd_me)
    if (ierr /= nf90_noerr) then
      message = trim(message)//"inquire_variable(mean_elev) failed: "//trim(nf90_strerror(ierr))
      return
    end if

    if (nd_af /= nd_me) then
      message = trim(message)//"area_frac and mean_elev have different rank; unsupported."
      ierr=20; return
    end if
    if (nd_af /= 2 .and. nd_af /= 3) then
      message = trim(message)//"expected rank-2 (spat2,band) or rank-3 (spat1,spat2,band)"
      ierr=20; return
    end if

    ! ---- read into canonical arrays ----
    
    !! Reads either:
    !!   rank-3 (spat1,spat2,band) into out3d(nx,ny,band)
    !!   rank-2 (spat2,band) into out3d(1,ny,band)  [requires nx_local==1]
    
    call read_elev_vars_to_canonical(ncid_eb, vid_af, vid_me, nd_af, info%space, info%snow%n_bands, &
                                     af, zmid, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! ---- compute weighted mean elevation + mask + AF sanity check ----
    do j = 1, info%space%ny_local
      do i = 1, info%space%nx_local
          
        domain%z_forcing(i,j) = sum(zmid(i,j,:) * af(i,j,:))

        ! if mean elevation first band is NA_VALUE, mask this grid cell
        domain%elev_mask(i,j) = zmid(i,j,1) == NA_VALUE_SP ! TODO check comparison against real
        if(.NOT.domain%elev_mask(i,j)) THEN ! only check area fraction sum to 1 if not NA_VALUE

          do ib=1,info%snow%n_bands
            domain%bands_info(i,j,ib)%num   = ib
            domain%bands_info(i,j,ib)%af    = af(i,j,ib)
            domain%bands_info(i,j,ib)%z_mid = zmid(i,j,ib)
          end do

          afsum = sum(af(i,j,:))
          if (abs(afsum - 1.0_wp) > 1.0e-2_wp) then
            write(message,'(a,2(i0,1x),a,f10.6)') trim(message)// &
              "AF sum != 1 at (i,j)= ", i, j, " sum=", afsum
            ierr=20; return
          end if

        endif  ! if masked

      end do
    end do

    ! --- close NetCDF file ---
    ierr = nf90_close(ncid_eb)
    if(ierr /= nf90_noerr) then
      message = trim(message)//"nf90_close failed: "//trim(nf90_strerror(ierr))// &
                " [file="//trim(eb_file)//"]"
      return
    endif

  end subroutine read_elevbands


  ! -------------------------------------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------------------------------------

  ! ----- private subroutine read_elev_vars_to_canonical (read elevation-band variables) ------------------------------

  subroutine read_elev_vars_to_canonical(ncid, vid_af, vid_me, ndims, space, n_bands, af, zmid, ierr, message)

  !----------------------------------------------------------------------------------------
  ! Read elevation-band variables (area_frac, mean_elev) from the elevation-band NetCDF
  ! file and map them into the model’s canonical (nx_local, ny_local, n_bands) layout.
  ! Supports both grid layout (spat1,spat2,band) and list layout (spat2,band → nx=1).
  ! Computes weighted mean forcing elevation for the local subdomain.
  !----------------------------------------------------------------------------------------

  implicit none

  integer(i4b), intent(in) :: ncid, vid_af, vid_me, ndims
  type(space_info), intent(in) :: space
  integer(i4b), intent(in) :: n_bands

  real(wp), intent(out) :: af(space%nx_local, space%ny_local, n_bands)
  real(wp), intent(out) :: zmid(space%nx_local, space%ny_local, n_bands)

  integer(i4b), intent(out) :: ierr
  character(*), intent(out) :: message

  integer(i4b) :: start3(3), count3(3)
  integer(i4b) :: start2(2), count2(2)
  real(wp)     :: tmp2_af(space%ny_local, n_bands)
  real(wp)     :: tmp2_me(space%ny_local, n_bands)

  ierr = 0
  message = "read_elev_vars_to_canonical/"

  if (ndims == 3) then

    start3 = (/ 1, space%y_start_global, 1 /)
    count3 = (/ space%nx_local, space%ny_local, n_bands /)

    ierr = nf90_get_var(ncid, vid_af, af, start=start3, count=count3)
    if (ierr /= nf90_noerr) then
      message = trim(message)//"nf90_get_var(area_frac) failed: "//trim(nf90_strerror(ierr))
      return
    endif

    ierr = nf90_get_var(ncid, vid_me, zmid, start=start3, count=count3)
    if (ierr /= nf90_noerr) then
      message = trim(message)//"nf90_get_var(mean_elev) failed: "//trim(nf90_strerror(ierr))
      return
    endif

  else

    ! rank-2 list case
    if (space%nx_local /= 1) then
      message = trim(message)//"rank-2 elevband file requires nx_local=1"
      ierr = 20; return
    endif

    start2 = (/ space%y_start_global, 1 /)
    count2 = (/ space%ny_local,       n_bands /)

    ierr = nf90_get_var(ncid, vid_af, tmp2_af, start=start2, count=count2)
    if (ierr /= nf90_noerr) then
      message = trim(message)//"nf90_get_var(area_frac list) failed: "//trim(nf90_strerror(ierr))
      return
    endif

    ierr = nf90_get_var(ncid, vid_me, tmp2_me, start=start2, count=count2)
    if (ierr /= nf90_noerr) then
      message = trim(message)//"nf90_get_var(mean_elev list) failed: "//trim(nf90_strerror(ierr))
      return
    endif

    af(1,:,:)   = tmp2_af
    zmid(1,:,:) = tmp2_me

  endif

  end subroutine

end module read_elevbands_module
