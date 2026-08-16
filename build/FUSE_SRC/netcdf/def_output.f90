MODULE DEF_OUTPUT_MODULE

  USE nrtype
  USE netcdf

  USE info_types,   only: fuse_info
  USE domain_types, only: coord_data

  USE fuse_globaldata, only: VAR_OBS, VAR_BAND, VAR_BASIN, VAR_REACH

  use iso_fortran_env, only: real32

  implicit none

  private
  public :: DEF_OUTPUT

contains

  SUBROUTINE DEF_OUTPUT(info, coords)

    USE metaoutput, only: VARDESCRIBE
    USE fuse_globaldata, only: FUSE_VERSION, FUSE_BUILDTIME, FUSE_GITBRANCH, FUSE_GITHASH
    USE metaoutput, only: NOUTVAR, VNAME, LNAME, VUNIT, VTYPE, isFlux
    USE model_defn, only: FNAME_NETCDF_RUNS
    USE fuse_fileManager, only: Q_ONLY
    USE multiforce, only: timeUnits
    USE fuse_globaldata, only: ncid_out

    implicit none

    type(fuse_info),    intent(in) :: info
    type(coord_data),   intent(in) :: coords

    ! locals
    integer(i4b) :: nPar, nObs, nSpat1, nSpat2, n_bands, n_seg
    integer(i4b) :: ierr, ivar, varid, varid_time, varid_lat, varid_lon, varid_band, varid_param, varid_seg
    integer(i4b) :: dim_time, dim_x, dim_y, dim_band, dim_par, dim_obs, dim_seg
    integer(i4b), dimension(3) :: dimids_basin
    integer(i4b), dimension(2) :: dimids_reach
    integer(i4b), dimension(4) :: dimids_band
    integer(i4b), dimension(4) :: dimids_par
    integer(i4b), dimension(2) :: dimids_obs

    real(real32), dimension(info%space%nx_local, &
                            info%space%ny_local) :: longitude
    real(real32), dimension(info%space%nx_local, &
                            info%space%ny_local) :: latitude
    
    integer(i4b), dimension(info%snow%n_bands)   :: band_i
    integer(i4b), dimension(info%config%nParam)  :: param_i
    integer(i4b) :: ib, ip

    integer(i4b) :: iOut

    real(real32), parameter                :: NA_VALUE_OUT = -9999._real32

    nPar    = info%config%nParam

    nObs    = info%space%nObs
    
    nSpat1  = info%space%nx_local  ! NOTE: local to rank (MPI parallelization)
    nSpat2  = info%space%ny_local
    
    n_seg   = info%space%n_seg
    n_bands = info%snow%n_bands

    call VARDESCRIBE()

    print *, 'Create NetCDF file for runs:'
    print *, trim(FNAME_NETCDF_RUNS)

    ! Create NetCDF-4 file (HDF5 container)
    ierr = nf90_create(trim(FNAME_NETCDF_RUNS), NF90_CLASSIC_MODEL, ncid_out)
    call handle_err(ierr)

    ! Dimensions
    ierr = nf90_def_dim(ncid_out, "time",  NF90_UNLIMITED, dim_time); call handle_err(ierr)
    ierr = nf90_def_dim(ncid_out, "band",  n_bands,        dim_band); call handle_err(ierr)
    ierr = nf90_def_dim(ncid_out, "param", nPar,           dim_par);  call handle_err(ierr)
    ierr = nf90_def_dim(ncid_out, "obs",   nObs,           dim_obs);  call handle_err(ierr)
    ierr = nf90_def_dim(ncid_out, "seg",   n_seg,          dim_seg);  call handle_err(ierr)
    ierr = nf90_def_dim(ncid_out, "x",     nSpat1,         dim_x);    call handle_err(ierr)
    ierr = nf90_def_dim(ncid_out, "y",     nSpat2,         dim_y);    call handle_err(ierr)

    dimids_basin = (/ dim_x, dim_y,           dim_time /)
    dimids_reach = (/               dim_seg,  dim_time /)
    dimids_band  = (/ dim_x, dim_y, dim_band, dim_time /)
    dimids_par   = (/ dim_x, dim_y, dim_par,  dim_time /)
    dimids_obs   = (/               dim_obs,  dim_time /)

    ! loop through desired output variables
    ! NOTE: Does not execute loop if nOutput=0
    do iOut = 1, info%config%nOutput

      ! find index of the variable name
      ivar = findloc(VNAME, info%config%outvar_names(iOut), dim=1)
      if (ivar == 0) cycle ! allow some names not in the metadata structure

      ! define variable
      select case( VTYPE(ivar) )
        case (VAR_OBS);   ierr = nf90_def_var(ncid_out, trim(VNAME(ivar)), NF90_FLOAT, dimids_obs,   varid)
        case (VAR_BAND);  ierr = nf90_def_var(ncid_out, trim(VNAME(ivar)), NF90_FLOAT, dimids_band,  varid)
        case (VAR_BASIN); ierr = nf90_def_var(ncid_out, trim(VNAME(ivar)), NF90_FLOAT, dimids_basin, varid)
        case (VAR_REACH); ierr = nf90_def_var(ncid_out, trim(VNAME(ivar)), NF90_FLOAT, dimids_reach, varid)
        case default; stop 'DEF_OUTPUT/Unable to identify variable type'
      end select
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

    ! Coordinate variables
    ierr = nf90_def_var(ncid_out, "time", NF90_FLOAT, (/dim_time/), varid_time);         call handle_err(ierr)
    ierr = nf90_put_att(ncid_out, varid_time, "units", trim(timeUnits));                 call handle_err(ierr)

    ierr = nf90_def_var(ncid_out, "latitude",  NF90_FLOAT, (/dim_x, dim_y/), varid_lat); call handle_err(ierr)
    ierr = nf90_put_att(ncid_out, varid_lat, "standard_name", "latitude");               call handle_err(ierr)
    ierr = nf90_put_att(ncid_out, varid_lat, "units", "degrees_north");                  call handle_err(ierr)
    
    ierr = nf90_def_var(ncid_out, "longitude", NF90_FLOAT, (/dim_x, dim_y/), varid_lon); call handle_err(ierr)
    ierr = nf90_put_att(ncid_out, varid_lon, "standard_name", "longitude");              call handle_err(ierr)
    ierr = nf90_put_att(ncid_out, varid_lon, "units", "degrees_east");                   call handle_err(ierr)

    ierr = nf90_def_var(ncid_out, "param", NF90_INT, (/dim_par/), varid_param);          call handle_err(ierr)
    ierr = nf90_put_att(ncid_out, varid_param, "units", "-");                            call handle_err(ierr)

    ierr = nf90_def_var(ncid_out, "band", NF90_INT, (/dim_band/), varid_band);           call handle_err(ierr)
    ierr = nf90_put_att(ncid_out, varid_band, "units", "-");                             call handle_err(ierr)

    ierr = nf90_def_var(ncid_out, "seg", NF90_INT, (/dim_seg/), varid_seg);              call handle_err(ierr)
    ierr = nf90_put_att(ncid_out, varid_seg, "units", "-");                              call handle_err(ierr)

    ! Global attributes
    ierr = nf90_put_att(ncid_out, NF90_GLOBAL, "software",        "FUSE");               call handle_err(ierr)
    ierr = nf90_put_att(ncid_out, NF90_GLOBAL, "fuse_version",    trim(FUSE_VERSION));   call handle_err(ierr)
    ierr = nf90_put_att(ncid_out, NF90_GLOBAL, "fuse_build_time", trim(FUSE_BUILDTIME)); call handle_err(ierr)
    ierr = nf90_put_att(ncid_out, NF90_GLOBAL, "fuse_git_branch", trim(FUSE_GITBRANCH)); call handle_err(ierr)
    ierr = nf90_put_att(ncid_out, NF90_GLOBAL, "fuse_git_hash",   trim(FUSE_GITHASH));   call handle_err(ierr)

    ! Leave define mode
    ierr = nf90_enddef(ncid_out); call handle_err(ierr)

    ! Write coordinate data
    latitude  = real(coords%lat_2d, kind(real32))
    longitude = real(coords%lon_2d, kind(real32))

    ierr = nf90_put_var(ncid_out, varid_lat, latitude);  call handle_err(ierr)
    ierr = nf90_put_var(ncid_out, varid_lon, longitude); call handle_err(ierr)

    band_i  = [(ib, ib=1,n_bands)]
    param_i = [(ip, ip=1,nPar)]

    ierr = nf90_put_var(ncid_out, varid_band,  band_i);  call handle_err(ierr)
    ierr = nf90_put_var(ncid_out, varid_param, param_i); call handle_err(ierr)
    
    ierr = nf90_put_var(ncid_out, varid_seg, coords%seg_id); call handle_err(ierr)

    ierr = nf90_close(ncid_out); call handle_err(ierr)

    print *, 'NetCDF file for model runs defined with dimensions', nSpat1, nSpat2, n_bands, nPar

  END SUBROUTINE DEF_OUTPUT

END MODULE DEF_OUTPUT_MODULE
