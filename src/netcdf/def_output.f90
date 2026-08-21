MODULE DEF_OUTPUT_MODULE

  USE nrtype
  USE netcdf

  USE info_types,   only: fuse_info
  USE domain_types, ONLY: domain_data

  USE fuse_globaldata, only: VAR_OBS, VAR_BAND, VAR_BASIN, VAR_REACH

  use iso_fortran_env, only: real32

  implicit none

  private
  public :: DEF_OUTPUT

contains

  SUBROUTINE DEF_OUTPUT(info, domain)

    USE fuse_globaldata,   only: FUSE_VERSION, FUSE_BUILDTIME, FUSE_GITBRANCH, FUSE_GITHASH
    USE fuse_globaldata,   only: do_mizuRoute
    USE fuse_globaldata,   only: ncid_out
    
    USE globaldata,        only: routeMethods ! mizuRoute
    
    USE handle_err_module, only: handle_err

    USE metaoutput,        only: VARDESCRIBE
    USE metaoutput,        only: NOUTVAR, VNAME, LNAME, VUNIT, VTYPE, isFlux

    USE init_mizuRoute,    only: route_method_name

    implicit none

    type(fuse_info),    intent(in) :: info
    type(domain_data),  intent(in) :: domain

    ! locals
    integer(i4b) :: nPar, nObs, nSpat1, nSpat2, n_bands, n_seg
    integer(i4b) :: ierr, ivar 
    integer(i4b) :: varid, varid_time, varid_time_bnds, varid_lat, varid_lon
    integer(i4b) :: varid_band, varid_param, varid_seg, varid_method
    integer(i4b) :: dim_time, dim_bnds, dim_x, dim_y 
    integer(i4b) :: dim_band, dim_par, dim_obs, dim_seg, dim_method
    integer(i4b), dimension(3) :: dimids_basin
    integer(i4b), dimension(3) :: dimids_reach
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

    character(len=32) :: attName
    integer(i4b)      :: iRoute

    character(len=32) :: subname
    subname="put_output/"

    nPar    = info%config%nParam

    nObs    = info%space%nObs
    
    nSpat1  = info%space%nx_local  ! NOTE: local to rank (MPI parallelization)
    nSpat2  = info%space%ny_local
    
    n_seg   = info%space%n_seg
    n_bands = info%snow%n_bands

    ! check the routing model is active
    do_mizuRoute = allocated(routeMethods)

    ! build metadata vectors for all variables
    call VARDESCRIBE() ! NOUTVAR, VNAME, LNAME, VUNIT, VTYPE, isFlux

    ! early return: no time-series output requested
    if ( .not. info%config%write_timeseries ) return

    print *, 'Create NetCDF file for runs:'
    print *, trim(info%files%fname_netcdf_runs)

    ! Create NetCDF-4 file (HDF5 container)
    ierr = nf90_create(trim(info%files%fname_netcdf_runs), NF90_CLASSIC_MODEL, ncid_out)
    call handle_err(ierr)

    ! Dimensions (land model)
    ierr = nf90_def_dim(ncid_out, "time",  NF90_UNLIMITED, dim_time); call handle_err(ierr)
    ierr = nf90_def_dim(ncid_out, "nbnds", 2,              dim_bnds); call handle_err(ierr)
    ierr = nf90_def_dim(ncid_out, "band",  n_bands,        dim_band); call handle_err(ierr)
    ierr = nf90_def_dim(ncid_out, "param", nPar,           dim_par);  call handle_err(ierr)
    ierr = nf90_def_dim(ncid_out, "obs",   nObs,           dim_obs);  call handle_err(ierr)
    ierr = nf90_def_dim(ncid_out, "x",     nSpat1,         dim_x);    call handle_err(ierr)
    ierr = nf90_def_dim(ncid_out, "y",     nSpat2,         dim_y);    call handle_err(ierr)

    dimids_basin = (/ dim_x, dim_y,           dim_time /)
    dimids_band  = (/ dim_x, dim_y, dim_band, dim_time /)
    dimids_par   = (/ dim_x, dim_y, dim_par,  dim_time /)
    dimids_obs   = (/               dim_obs,  dim_time /)
    
    ! Dimensions (routing model)
    if ( do_mizuRoute ) then
    
     ierr = nf90_def_dim(ncid_out, "seg",    n_seg,              dim_seg);     call handle_err(ierr)
     ierr = nf90_def_dim(ncid_out, "method", size(routeMethods), dim_method);  call handle_err(ierr)
    
     dimids_reach = (/dim_method, dim_seg, dim_time /)

    end if

    ! loop through desired output variables
    ! NOTE: Does not execute loop if nOutput=0
    do iOut = 1, info%config%nOutput

      ! find index of the variable name
      ivar = findloc(VNAME, info%config%outvar_names(iOut), dim=1)
      if (ivar == 0) cycle ! allow some names not in the metadata structure

      if ( VTYPE(ivar) == VAR_REACH .and. .not. do_mizuRoute)  &
      stop trim(subname)//": output variable "//trim(VNAME(ivar))//" requires mizuRoute"

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
    ierr = nf90_def_var(ncid_out, "time", NF90_DOUBLE, (/dim_time/), varid_time);        call handle_err(ierr)
    ierr = nf90_put_att(ncid_out, varid_time, "axis", "T");                              call handle_err(ierr, trim(subname))
    ierr = nf90_put_att(ncid_out, varid_time, "bounds", "time_bnds");                    call handle_err(ierr, trim(subname))
    ierr = nf90_put_att(ncid_out, varid_time, "standard_name", "time");                  call handle_err(ierr, trim(subname))
    ierr = nf90_put_att(ncid_out, varid_time, "units", trim(info%time%units));           call handle_err(ierr, trim(subname))
    ierr = nf90_put_att(ncid_out, varid_time, "long_name", "midpoint of time interval"); call handle_err(ierr, trim(subname))

    ierr = nf90_def_var(ncid_out, "time_bnds", NF90_DOUBLE, (/dim_bnds, dim_time/), varid_time_bnds); call handle_err(ierr)
    ierr = nf90_put_att(ncid_out, varid_time_bnds, "units", trim(info%time%units));      call handle_err(ierr)

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

    if (do_mizuRoute) then

      ierr = nf90_def_var(ncid_out, "seg", NF90_INT, (/dim_seg/), varid_seg);            call handle_err(ierr)
      ierr = nf90_put_att(ncid_out, varid_seg, "units", "-");                            call handle_err(ierr)

      ierr = nf90_def_var(ncid_out, "method", NF90_INT, (/dim_method/), varid_method);   call handle_err(ierr)
      ierr = nf90_put_att(ncid_out, varid_method, "long_name", "routing method");        call handle_err(ierr)
      ierr = nf90_put_att(ncid_out, varid_method, "source", "mizuRoute");                call handle_err(ierr)

      ! add decription of routing methods
      do iRoute = 1, size(routeMethods)
        write(attName,'("method_",I0)') routeMethods(iRoute)
        ierr = nf90_put_att(ncid_out, varid_method, trim(attName), &
                            route_method_name(routeMethods(iRoute)))
        call handle_err(ierr, trim(subname)//":nf90_put_att(method:"//trim(attName)//")")                 
      end do

    end if

    ! Global attributes
    ierr = nf90_put_att(ncid_out, NF90_GLOBAL, "software",        "FUSE");               call handle_err(ierr)
    ierr = nf90_put_att(ncid_out, NF90_GLOBAL, "fuse_version",    trim(FUSE_VERSION));   call handle_err(ierr)
    ierr = nf90_put_att(ncid_out, NF90_GLOBAL, "fuse_build_time", trim(FUSE_BUILDTIME)); call handle_err(ierr)
    ierr = nf90_put_att(ncid_out, NF90_GLOBAL, "fuse_git_branch", trim(FUSE_GITBRANCH)); call handle_err(ierr)
    ierr = nf90_put_att(ncid_out, NF90_GLOBAL, "fuse_git_hash",   trim(FUSE_GITHASH));   call handle_err(ierr)

    ! Leave define mode
    ierr = nf90_enddef(ncid_out); call handle_err(ierr)

    ! Write coordinate data
    latitude  = real(domain%coords%lat_2d, kind(real32))
    longitude = real(domain%coords%lon_2d, kind(real32))

    ierr = nf90_put_var(ncid_out, varid_lat, latitude);  call handle_err(ierr)
    ierr = nf90_put_var(ncid_out, varid_lon, longitude); call handle_err(ierr)

    band_i  = [(ib, ib=1,n_bands)]
    param_i = [(ip, ip=1,nPar)]

    ierr = nf90_put_var(ncid_out, varid_band,  band_i);  call handle_err(ierr)
    ierr = nf90_put_var(ncid_out, varid_param, param_i); call handle_err(ierr)
    
    if (do_mizuRoute) then
      ierr = nf90_put_var(ncid_out, varid_seg,    domain%reach%seg_id); call handle_err(ierr)
      ierr = nf90_put_var(ncid_out, varid_method, routeMethods);        call handle_err(ierr)
    endif

    ierr = nf90_close(ncid_out); call handle_err(ierr)

    print *, 'NetCDF file for model runs defined with dimensions', nSpat1, nSpat2, n_bands, nPar

  END SUBROUTINE DEF_OUTPUT

END MODULE DEF_OUTPUT_MODULE
