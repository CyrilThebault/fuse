module manage_legacy_data

  use nrtype, only: i4b

  implicit none

  private

  public :: deallocate_legacy_data

contains

subroutine deallocate_legacy_data(ierr, message)

  use multiforce, only: aValid

  implicit none

  integer(i4b), intent(out) :: ierr
  character(*), intent(out) :: message

  ierr    = 0
  message = 'deallocate_legacy_data/'

  DEALLOCATE(aValid, stat=ierr)
  if(ierr/=0)then
    message=trim(message)//'unable to deallocate space for catchment modeling'
    return
  endif

end subroutine deallocate_legacy_data

end module manage_legacy_data
