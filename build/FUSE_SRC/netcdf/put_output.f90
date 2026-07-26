module put_output_module

  use nrtype
  use work_types, only: fuse_work
  use iso_fortran_env, only: real32

  use netcdf, only: &
      NF90_WRITE, NF90_NOERR, &
      nf90_open, nf90_close, nf90_inq_varid, nf90_put_var

  implicit none
  private
  public :: put_output

contains

  subroutine put_output(fuseStruct, istart_sim, istart_in, numtim)

  ! -------------------------------------------------------------------------------------
  ! Creator:
  ! --------
  ! Nans Addor, based on Martyn Clark's 2007 PUT_OUTPUT
  ! Modified by Martyn Clark to use the elevation band dimension and add parameter derivatives, 12/2025
  ! Modified by Martyn Clark to use output buffers in fuseStruct
  ! -------------------------------------------------------------------------------------
  ! Purpose:
  ! --------
  !   Write a 3D (or 4D) data structure to the NetCDF output file (chunk write)
  ! -------------------------------------------------------------------------------------

  ! subroutines
  use varextract_module, only: varextract_3d

  ! metadata / config
  use model_defn,    only: fname_netcdf_runs
  use metaoutput,    only: noutvar, vname, isband
  use multiparam,    only: numpar
  use multibands,    only: mbands_var_4d, n_bands
  use multiforce,    only: time_steps, nspat1, nspat2
  use fuse_filemanager, only: q_only

  ! global
  use globaldata,    only: ncid_out

  implicit none

  ! input
  type(fuse_work), intent(in) :: fuseStruct
  integer(i4b),    intent(in) :: istart_sim
  integer(i4b),    intent(in) :: istart_in
  integer(i4b),    intent(in) :: numtim

  ! locals
  logical(lgt) :: write_var
  integer(i4b) :: ierr
  integer(i4b) :: ivar
  integer(i4b) :: ivar_id

  integer(i4b), dimension(3) :: start3, count3
  integer(i4b), dimension(4) :: start4_band, count4_band
  integer(i4b), dimension(4) :: start4_param, count4_param

  real(real32), dimension(nspat1, nspat2, numtim)              :: avar_3d

  real(real32), dimension(nspat1, nspat2, n_bands, numtim)     :: avar_4d_band
  ! placeholder for future param-derivative write
  real(real32), dimension(nspat1, nspat2, numpar,  numtim)     :: avar_4d_param

  real(real32), dimension(numtim) :: time_steps_sub

  character(len=32) :: subname
  subname="put_output.f90"

  ! -----------------------------------------------------------------------------
  ! dimension lists (Fortran nf90 uses 1-based indices)
  start3 = (/1, 1, istart_sim/)
  count3 = (/nspat1, nspat2, numtim/)

  start4_band = (/1, 1, 1, istart_sim/)
  count4_band = (/nspat1, nspat2, n_bands, numtim/)

  start4_param = (/1, 1, 1, istart_sim/)
  count4_param = (/nspat1, nspat2, numpar,  numtim/)

  ! open file (already defined elsewhere via DEF_OUTPUT)
  ierr = nf90_open(trim(fname_netcdf_runs), NF90_WRITE, ncid_out)
  call handle_err(ierr, trim(subname)//":nf90_open")

  ! loop through variables with time-varying model output
  do ivar = 1, noutvar

    ! optional "Q_ONLY" filter
    if (q_only) then
      select case (trim(vname(ivar)))
        case ('q_instnt', 'q_routed');  write_var = .true.
        case default;                   write_var = .false.
      end select
    end if

    if (.not. write_var) cycle

    ! get var id
    ierr = nf90_inq_varid(ncid_out, trim(vname(ivar)), ivar_id)
    call handle_err(ierr, trim(subname)//":nf90_inq_varid:"//trim(vname(ivar)))

    if (.not. isband(ivar)) then

      ! 3-d variable -- extract from the output buffers in fuseStruct%chunk
      call varextract_3d(fuseStruct%chunk, vname(ivar), nspat1, nspat2, numtim, avar_3d)

      ierr = nf90_put_var(ncid_out, ivar_id, avar_3d, start=start3, count=count3)
      call handle_err(ierr, trim(subname)//":nf90_put_var(3d):"//trim(vname(ivar)))

    else

      ! 4-d elevation band variable (stored in MBANDS_VAR_4d)
      select case (trim(vname(ivar)))
        case ('swe_z');     avar_4d_band = mbands_var_4d(:,:,:,1:numtim)%swe
        case ('snwacml_z'); avar_4d_band = mbands_var_4d(:,:,:,1:numtim)%snowaccmltn
        case ('snwmelt_z'); avar_4d_band = mbands_var_4d(:,:,:,1:numtim)%snowmelt
        case default;       stop trim(subname)//":unknown band var:"//trim(vname(ivar))
      end select

      ierr = nf90_put_var(ncid_out, ivar_id, avar_4d_band, start=start4_band, count=count4_band)
      call handle_err(ierr, trim(subname)//":nf90_put_var(4d band):"//trim(vname(ivar)))

    end if

    ! future: param-derivative writes would go here using count4_param/start4_param
    ! e.g. name = trim(vname(ivar))//'__dFlux_dParam'

  end do

  ! write time
  time_steps_sub = real(time_steps(istart_in:(istart_in + numtim - 1)), kind(real32))

  ierr = nf90_inq_varid(ncid_out, 'time', ivar_id)
  call handle_err(ierr, trim(subname)//":nf90_inq_varid:time")

  ierr = nf90_put_var(ncid_out, ivar_id, time_steps_sub, start=(/istart_sim/), count=(/numtim/))
  call handle_err(ierr, trim(subname)//":nf90_put_var:time")

  ! close
  ierr = nf90_close(ncid_out)
  call handle_err(ierr, trim(subname)//":nf90_close")

  end subroutine put_output

end module put_output_module
