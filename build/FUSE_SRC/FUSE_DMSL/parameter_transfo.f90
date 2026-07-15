! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Cyril Thébault, 2026
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Apply parameter transformations for optimization.
!
! Parameters remain in physical space inside FUSE. Optimizers may operate in a
! transformed search space. This module provides conversions between both spaces.
! ---------------------------------------------------------------------------------------

module parameter_transform_module

  use nrtype, only : MSP, I4B

  implicit none
  private

  ! Transformation codes used in the PARVTN column of zConstraints.
  integer(I4B), parameter, public :: TRANS_NONE  = 0
  integer(I4B), parameter, public :: TRANS_LOG10 = 1
  integer(I4B), parameter, public :: TRANS_LN    = 2

  public :: to_search_space
  public :: to_physical_space
  public :: vector_to_search_space
  public :: vector_to_physical_space
  public :: validate_transform

contains

  ! -------------------------------------------------------------------------------------
  ! Transform one physical parameter to the optimizer search space.
  ! -------------------------------------------------------------------------------------
  pure function to_search_space(x, transform_code) result(z)

    real(MSP), intent(in)    :: x
    integer(I4B), intent(in) :: transform_code
    real(MSP)                :: z

    select case (transform_code)

    case (TRANS_NONE)
      z = x

    case (TRANS_LOG10)
      z = log10(x)

    case (TRANS_LN)
      z = log(x)

    case default
      ! Unknown codes must be rejected by validate_transform before this function is used.
      z = x

    end select

  end function to_search_space


  ! -------------------------------------------------------------------------------------
  ! Transform one optimizer value back to the physical parameter space.
  ! -------------------------------------------------------------------------------------
  pure function to_physical_space(z, transform_code) result(x)

    real(MSP), intent(in)    :: z
    integer(I4B), intent(in) :: transform_code
    real(MSP)                :: x

    select case (transform_code)

    case (TRANS_NONE)
      x = z

    case (TRANS_LOG10)
      x = 10.0_MSP**z

    case (TRANS_LN)
      x = exp(z)

    case default
      ! Unknown codes must be rejected by validate_transform before this function is used.
      x = z

    end select

  end function to_physical_space


  ! -------------------------------------------------------------------------------------
  ! Transform a complete parameter vector to the optimizer search space.
  ! -------------------------------------------------------------------------------------
  subroutine vector_to_search_space(physical_values, transform_codes, search_values, &
                                    ierr, message)

    real(MSP), intent(in)                 :: physical_values(:)
    integer(I4B), intent(in)              :: transform_codes(:)
    real(MSP), intent(out)                :: search_values(:)
    integer(I4B), intent(out)             :: ierr
    character(len=*), intent(out)         :: message

    integer(I4B) :: i
    integer(I4B) :: npar

    ierr    = 0
    message = ''

    npar = size(physical_values)

    if (size(transform_codes) /= npar .or. size(search_values) /= npar) then
      ierr = 10
      message = 'Inconsistent array sizes in vector_to_search_space'
      return
    end if

    do i = 1, npar

      select case (transform_codes(i))

      case (TRANS_NONE)
        continue

      case (TRANS_LOG10, TRANS_LN)

        if (physical_values(i) <= 0.0_MSP) then
          ierr = 11
          write(message,'(a,i0,a)') &
            'Cannot apply logarithmic transformation to parameter ', i, &
            ': value must be strictly positive.'
          return
        end if

      case default
        ierr = 12
        write(message,'(a,i0,a,i0)') &
          'Unknown transformation code ', transform_codes(i), &
          ' for parameter ', i
        return

      end select

      search_values(i) = to_search_space(physical_values(i), transform_codes(i))

    end do

  end subroutine vector_to_search_space


  ! -------------------------------------------------------------------------------------
  ! Transform a complete optimizer vector back to physical parameter space.
  ! -------------------------------------------------------------------------------------
  subroutine vector_to_physical_space(search_values, transform_codes, physical_values, &
                                      ierr, message)

    real(MSP), intent(in)                 :: search_values(:)
    integer(I4B), intent(in)              :: transform_codes(:)
    real(MSP), intent(out)                :: physical_values(:)
    integer(I4B), intent(out)             :: ierr
    character(len=*), intent(out)         :: message

    integer(I4B) :: i
    integer(I4B) :: npar

    ierr    = 0
    message = ''

    npar = size(search_values)

    if (size(transform_codes) /= npar .or. size(physical_values) /= npar) then
      ierr = 20
      message = 'Inconsistent array sizes in vector_to_physical_space'
      return
    end if

    do i = 1, npar

      select case (transform_codes(i))

      case (TRANS_NONE, TRANS_LOG10, TRANS_LN)
        physical_values(i) = &
          to_physical_space(search_values(i), transform_codes(i))

      case default
        ierr = 21
        write(message,'(a,i0,a,i0)') &
          'Unknown transformation code ', transform_codes(i), &
          ' for parameter ', i
        return

      end select

    end do

  end subroutine vector_to_physical_space


  ! -------------------------------------------------------------------------------------
  ! Validate a transformation and its physical bounds.
  ! -------------------------------------------------------------------------------------
  subroutine validate_transform(parname, lower, default_value, upper, &
                                transform_code, ierr, message)

    character(len=*), intent(in)          :: parname
    real(MSP), intent(in)                 :: lower
    real(MSP), intent(in)                 :: default_value
    real(MSP), intent(in)                 :: upper
    integer(I4B), intent(in)              :: transform_code
    integer(I4B), intent(out)             :: ierr
    character(len=*), intent(out)         :: message

    ierr    = 0
    message = ''

    if (lower > upper) then
      ierr = 1
      write(message,'(a,a)') &
        'Invalid parameter bounds for ', trim(parname)
      return
    end if

    if (default_value < lower .or. default_value > upper) then
      ierr = 2
      write(message,'(a,a)') &
        'Default value outside bounds for ', trim(parname)
      return
    end if

    select case (transform_code)

    case (TRANS_NONE)
      continue

    case (TRANS_LOG10, TRANS_LN)

      if (lower <= 0.0_MSP .or. default_value <= 0.0_MSP .or. &
          upper <= 0.0_MSP) then

        ierr = 3
        write(message,'(a,a,a)') &
          'Logarithmic transformation requires strictly positive values for ', &
          trim(parname), '.'
        return

      end if

    case default

      ierr = 4
      write(message,'(a,i0,a,a)') &
        'Unknown transformation code ', transform_code, &
        ' for parameter ', trim(parname)
      return

    end select

  end subroutine validate_transform

end module parameter_transform_module
