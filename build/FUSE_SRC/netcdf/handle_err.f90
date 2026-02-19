subroutine handle_err(ierr, where)
  use nrtype,  only: i4b
  use netcdf,  only: NF90_NOERR, nf90_strerror
  implicit none

  integer(i4b), intent(in) :: ierr
  character(len=*), intent(in), optional :: where

  if (ierr /= NF90_NOERR) then
    if (present(where)) then
      write(*,'(a,1x,a)') 'NetCDF error in '//trim(where)//':', trim(nf90_strerror(ierr))
    else
      write(*,'(a)') trim(nf90_strerror(ierr))
    end if
    stop 1
  end if
end subroutine handle_err
