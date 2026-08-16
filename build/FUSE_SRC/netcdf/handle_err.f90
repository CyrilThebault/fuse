module handle_err_MODULE

  use nrtype,  only: i4b
  use netcdf,  only: NF90_NOERR, nf90_strerror

  implicit none

  private
  public :: handle_err

contains

  subroutine handle_err(ierr, string)
  implicit none

  integer(i4b), intent(in) :: ierr
  character(len=*), intent(in), optional :: string

    if (ierr /= NF90_NOERR) then
      if (present(string)) then
        write(*,'(a,1x,a)') 'NetCDF error in '//trim(string)//':', trim(nf90_strerror(ierr))
      else
        write(*,'(a)') trim(nf90_strerror(ierr))
      end if
      stop 1
    end if
  end subroutine handle_err

end module handle_err_MODULE
