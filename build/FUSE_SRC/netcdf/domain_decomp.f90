module domain_decomp_module

  use info_types, only: fuse_info

  implicit none

  private
  public :: get_domain_decomp_indices

contains

  ! -------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------

  ! ----- get indices to decompose the spatial domain -----------------------------------
  ! 1) Determine global run mode (grid vs catchment)
  ! 2) Apply MPI decomposition (y dimension) and store local dims + offsets

  subroutine get_domain_decomp_indices(info)
    implicit none
    type(fuse_info), intent(inout) :: info
   
    associate(&
        nx_global      => info%space%nx_global,      &
        ny_global      => info%space%ny_global,      &
        nx_local       => info%space%nx_local,       &
        ny_local       => info%space%ny_local,       &
        y_start_global => info%space%y_start_global, &
        y_end_global   => info%space%y_end_global,   &
        mpi_enabled    => info%mpi%enabled,          &
        nproc          => info%mpi%nproc,            &
        rank           => info%mpi%rank   )
   
    ! Copy globals
    nx_local = nx_global
    ny_local = ny_global
    y_start_global = 1
   
    ! Get indices for split dimensions
    if(mpi_enabled .and. nproc>1) then
      call split_1d(ny_global, rank, nproc, &  ! input
                    y_start_global, ny_local)  ! output
    endif
    y_end_global = y_start_global + ny_local - 1
   
    end associate
  end subroutine get_domain_decomp_indices

  ! -------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------

  ! ----- split the dimensions for each MPI rank ----------------------------------------
  ! Purpose: Split domain to allow for MPI.
  !          Given rank, nproc, and n_global, provide start and n_local indices
  ! Creator: Ethan Gutmann, 2020
  ! Modified by Martyn Clark to simplify code and input/output, 12/2025

  subroutine split_1d(n_global, rank, nproc, start, n_local, verbose)
  use nrtype
  implicit none
  integer(i4b), intent(in)  :: n_global, rank, nproc
  logical(lgt), intent(in), optional :: verbose
  integer(i4b), intent(out) :: start, n_local

  integer(i4b) :: base, extra
  logical(lgt) :: talk

  talk = .false.; if(present(verbose)) talk = verbose

  ! --- sanity checks ---
  if(nproc   <= 0)                    stop "split_1d: nproc must be > 0"
  if(rank     < 0 .or. rank >= nproc) stop "split_1d: rank out of range"
  if(n_global < 1)                    stop "split_1d: n_global must be >= 1"

  base  = n_global / nproc                      ! floor(n_global / nproc) rows per rank
  extra = mod(n_global, nproc)                  ! remainder; first 'extra' ranks get +1 row

  n_local = base + merge(1, 0, rank < extra)    ! add 1 row for ranks 0..extra-1
  start   = rank*base + min(rank, extra) + 1    ! shift start by #extra rows assigned before this rank

  if(talk) then
    write(*,'(a,i0,a,i0)') "split_1d: nproc=", nproc, " rank  =", rank
    write(*,'(a,i0,a,i0)') "split_1d: base =", base,  " extra =", extra
    write(*,'(a,i0,a,i0)') "split_1d: start=", start, " nLocal=", n_local
    write(*,'(a,i0,a,i0)') "split_1d: global rows=", start, ":", start+n_local-1
  endif
  end subroutine split_1d

end module domain_decomp_module
