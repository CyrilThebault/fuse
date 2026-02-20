MODULE DEF_OUTPUT_MODULE
  USE nrtype
  USE netcdf
  implicit none
  private
  public :: DEF_OUTPUT

contains

  SUBROUTINE DEF_OUTPUT(nSpat1,nSpat2,n_bands,NUMPAR)

    USE metaoutput, only: VARDESCRIBE
    USE globaldata, only: FUSE_VERSION, FUSE_BUILDTIME, FUSE_GITBRANCH, FUSE_GITHASH
    USE metaoutput, only: NOUTVAR, VNAME, LNAME, VUNIT, isBand, isFlux
    USE model_defn, only: FNAME_NETCDF_RUNS
    USE fuse_fileManager, only: Q_ONLY
    USE multiforce, only: latitude,longitude, timeUnits
    USE globaldata, only: ncid_out

    implicit none

    integer(i4b), intent(in) :: nSpat1, nSpat2, n_bands, NUMPAR

    ! locals
    integer(i4b) :: ierr, ivar, varid, varid_time, varid_lat, varid_lon, varid_band, varid_param
    integer(i4b) :: dim_time, dim_lon, dim_lat, dim_band, dim_par
    integer(i4b), dimension(3) :: dimids_3
    integer(i4b), dimension(4) :: dimids_band
    integer(i4b), dimension(4) :: dimids_par

    logical(lgt) :: write_var

    real(msp), dimension(nspat1) :: longitude_msp
    real(msp), dimension(nspat2) :: latitude_msp
    real(msp), parameter         :: NA_VALUE_OUT = -9999._msp

    integer(i4b), dimension(n_bands) :: band_i
    integer(i4b), dimension(NUMPAR)  :: param_i
    integer(i4b) :: ib, ip

    call VARDESCRIBE()

    print *, 'Create NetCDF file for runs:'
    print *, trim(FNAME_NETCDF_RUNS)

    ! Create NetCDF-4 file (HDF5 container)
    ierr = nf90_create(trim(FNAME_NETCDF_RUNS), NF90_CLASSIC_MODEL, ncid_out)
    call handle_err(ierr)

    ! Dimensions
    ierr = nf90_def_dim(ncid_out, "time", NF90_UNLIMITED, dim_time); call handle_err(ierr)
    ierr = nf90_def_dim(ncid_out, "band", n_bands, dim_band);        call handle_err(ierr)
    ierr = nf90_def_dim(ncid_out, "param", NUMPAR, dim_par);         call handle_err(ierr)
    ierr = nf90_def_dim(ncid_out, "longitude", nSpat1, dim_lon);     call handle_err(ierr)
    ierr = nf90_def_dim(ncid_out, "latitude",  nSpat2, dim_lat);     call handle_err(ierr)

    dimids_3    = (/ dim_lon, dim_lat, dim_time /)
    dimids_band = (/ dim_lon, dim_lat, dim_band, dim_time /)
    dimids_par  = (/ dim_lon, dim_lat, dim_par,  dim_time /)

    ! Time-varying output vars
    do ivar = 1, NOUTVAR

      if (Q_ONLY) then
        write_var = .false.
        if (trim(VNAME(ivar)) == "q_instnt") write_var = .true.
        if (trim(VNAME(ivar)) == "q_routed") write_var = .true.
        if (.not. write_var) cycle
      end if

      if (isBand(ivar)) then
        ierr = nf90_def_var(ncid_out, trim(VNAME(ivar)), NF90_FLOAT, dimids_band, varid)
      else
        ierr = nf90_def_var(ncid_out, trim(VNAME(ivar)), NF90_FLOAT, dimids_3, varid)
      end if
      call handle_err(ierr)
 
      ! Attributes
      ierr = nf90_put_att(ncid_out, varid, "long_name", trim(LNAME(ivar))); call handle_err(ierr)
      ierr = nf90_put_att(ncid_out, varid, "units",     trim(VUNIT(ivar))); call handle_err(ierr)
      ierr = nf90_put_att(ncid_out, varid, "_FillValue", NA_VALUE_OUT);     call handle_err(ierr)

      ! Optional: parameter sensitivity var for each flux
      if (isFlux(ivar)) then
        ierr = nf90_def_var(ncid_out, trim(VNAME(ivar))//"__dFlux_dParam", NF90_FLOAT, dimids_par, varid)
        call handle_err(ierr)
        ierr = nf90_put_att(ncid_out, varid, "_FillValue", NA_VALUE_OUT); call handle_err(ierr)
      end if

    end do  ! looping through variables

    print*, 'coordinate variables'

    ! Coordinate variables
    ierr = nf90_def_var(ncid_out, "time", NF90_FLOAT, (/dim_time/), varid_time); call handle_err(ierr)
    ierr = nf90_put_att(ncid_out, varid_time, "units", trim(timeUnits));         call handle_err(ierr)

    ierr = nf90_def_var(ncid_out, "latitude", NF90_FLOAT, (/dim_lat/), varid_lat); call handle_err(ierr)
    ierr = nf90_put_att(ncid_out, varid_lat, "units", "degreesN");                 call handle_err(ierr)
    ierr = nf90_put_att(ncid_out, varid_lat, "axis",  "Y");                        call handle_err(ierr)

    ierr = nf90_def_var(ncid_out, "longitude", NF90_FLOAT, (/dim_lon/), varid_lon); call handle_err(ierr)
    ierr = nf90_put_att(ncid_out, varid_lon, "units", "degreesE");                  call handle_err(ierr)
    ierr = nf90_put_att(ncid_out, varid_lon, "axis",  "X");                         call handle_err(ierr)

    ierr = nf90_def_var(ncid_out, "param", NF90_INT, (/dim_par/), varid_param); call handle_err(ierr)
    ierr = nf90_put_att(ncid_out, varid_param, "units", "-");                   call handle_err(ierr)

    ierr = nf90_def_var(ncid_out, "band", NF90_INT, (/dim_band/), varid_band);  call handle_err(ierr)
    ierr = nf90_put_att(ncid_out, varid_band, "units", "-");                    call handle_err(ierr)

    ! Global attributes
    ierr = nf90_put_att(ncid_out, NF90_GLOBAL, "software",        "FUSE");               call handle_err(ierr)
    ierr = nf90_put_att(ncid_out, NF90_GLOBAL, "fuse_version",    trim(FUSE_VERSION));   call handle_err(ierr)
    ierr = nf90_put_att(ncid_out, NF90_GLOBAL, "fuse_build_time", trim(FUSE_BUILDTIME)); call handle_err(ierr)
    ierr = nf90_put_att(ncid_out, NF90_GLOBAL, "fuse_git_branch", trim(FUSE_GITBRANCH)); call handle_err(ierr)
    ierr = nf90_put_att(ncid_out, NF90_GLOBAL, "fuse_git_hash",   trim(FUSE_GITHASH));   call handle_err(ierr)

    ! Leave define mode
    ierr = nf90_enddef(ncid_out); call handle_err(ierr)

    print*, 'coordinate data'

    print*, nSpat1, nSpat2
    print*, latitude
    print*, longitude

    ! Write coordinate data
    latitude_msp  = latitude
    longitude_msp = longitude

    print*, 'hello1'

    ierr = nf90_put_var(ncid_out, varid_lat, latitude_msp);  call handle_err(ierr)
    ierr = nf90_put_var(ncid_out, varid_lon, longitude_msp); call handle_err(ierr)

    print*, 'hello1'
    band_i  = [(ib, ib=1,n_bands)]
    param_i = [(ip, ip=1,NUMPAR)]

    print*, 'hello1'
    ierr = nf90_put_var(ncid_out, varid_band,  band_i);  call handle_err(ierr)
    ierr = nf90_put_var(ncid_out, varid_param, param_i); call handle_err(ierr)
    
    ierr = nf90_close(ncid_out); call handle_err(ierr)

    print *, 'NetCDF file for model runs defined with dimensions', nSpat1, nSpat2, n_bands, NUMPAR

  END SUBROUTINE DEF_OUTPUT

END MODULE DEF_OUTPUT_MODULE
