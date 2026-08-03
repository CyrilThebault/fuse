module sce_callback_context
  use nrtype,      only: I4B
  use info_types,  only: fuse_info
  use data_types,  only: domain_data
  use work_types,  only: fuse_work
  implicit none
  private
  public :: ctx
  public :: set_sce_context, clear_sce_context

  type :: sce_context
     type(fuse_info),   pointer :: info   => null()
     type(fuse_work),   pointer :: work   => null()
     type(domain_data), pointer :: domain => null()
     integer(I4B),      allocatable :: transform_codes(:)
  end type sce_context

  type(sce_context), save :: ctx

contains

  subroutine set_sce_context(info, work, domain, transform_codes)
    type(fuse_info),   target, intent(inout) :: info
    type(fuse_work),   target, intent(inout) :: work
    type(domain_data), target, intent(inout) :: domain
    integer(I4B),              intent(in)    :: transform_codes(:)
    ctx%info     => info
    ctx%work     => work
    ctx%domain   => domain

    if (allocated(ctx%transform_codes)) deallocate(ctx%transform_codes)

    allocate(ctx%transform_codes(size(transform_codes)))
    ctx%transform_codes = transform_codes

  end subroutine

  subroutine clear_sce_context()
    nullify(ctx%info, ctx%work, ctx%domain)
    if (allocated(ctx%transform_codes)) then
      deallocate(ctx%transform_codes)
    end if
  end subroutine

end module sce_callback_context
