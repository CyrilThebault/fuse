PROGRAM DISTRIBUTED_DRIVER
! ---------------------------------------------------------------------------------------
! Creators:
! Martyn Clark, 2011
! Modified by Brian Henn to include snow model, 6/2013
! Modified by Nans Addor to include distributed modeling, 9/2016
! Modified by Nans Addor to re-enable catchment-scale modeling, 4/2017
! Modified by Martyn Clark to modularize and simplify CLI, 12/2025
! ---------------------------------------------------------------------------------------
! Purpose:
! Driver program to run FUSE with a snow module as either at the catchment-scale or
! at the grid-scale
! ---------------------------------------------------------------------------------------
! data types
USE nrtype                                                ! variable types, etc.
USE info_types, only: cli_options                         ! command line interface options
USE info_types, only: fuse_info                           ! info structure (includes "everything")
USE work_types, only: fuse_work                           ! structures that depend on nState/nPar 
USE domain_types, only: domain_data                       ! domain data

! data
USE fuse_globaldata, only: isPrint
USE fuse_globaldata, only: ncid_out

! model setup: external subroutines/functions
USE netcdf                                                       ! NetCDF library
USE parse_command_args_MODULE, only: parse_command_args          ! parse command line arguments
USE setup_domain_module, only: setup_domain                      ! initialize the model domain
USE setup_model_definition_module, only: setup_model_definition  ! setup the FUSE model configuration
USE alloc_scratch_module, only: init_fuse_work                   ! initialze work structure

! model run: external subroutines/functions
USE get_fparam_module, only: GET_PRE_PARAM, GET_SCE_PARAM ! read parameters from netcdf file
USE sce_driver_MODULE, only: sce_driver                   ! SCE optimization

! model simulation modules
USE fuse_evaluate_module, only: fuse_evaluate             ! run model and compute performance metrics

#ifdef __MPI__
  use mpi
#endif

IMPLICIT NONE

! error control
integer(i4b)                           :: err              ! error code
character(len=1024)                    :: message          ! error message

! command line arguments
type(cli_options)                      :: cli_opts         ! command line argument options

! parameter set; parameter bounds
REAL(WP), DIMENSION(:), ALLOCATABLE    :: BL      ! vector of lower parameter bounds
REAL(WP), DIMENSION(:), ALLOCATABLE    :: BU      ! vector of upper parameter bounds
REAL(WP), DIMENSION(:), ALLOCATABLE    :: APAR    ! model parameter set

! function  evaluation
REAL(WP)                               :: METRIC_VAL      ! sim-obs differences

! model output
LOGICAL(LGT)                           :: OUTPUT_FLAG     ! .TRUE. = write time series output
INTEGER(I4B)                           :: ONEMOD=1        ! just specify one model

! global domain data
type(fuse_info)                        :: info            ! includes "everything"
type(fuse_work)                        :: work            ! structures that depend on nState/nPar
type(domain_data)                      :: domain          ! 3d/4d output buffers

! ---------------------------------------------------------------------------------------
! ----- model preliminaries (initialize) ------------------------------------------------
! ---------------------------------------------------------------------------------------

! ----- initialize MPI ------------------------------------------------------------------

#ifdef __MPI__
  info%mpi%enabled = .true.
  call MPI_Init(err); call MPI_check(err, "MPI_Init")
  call MPI_Comm_rank(MPI_COMM_WORLD, info%mpi%rank, err);  call MPI_check(err, "MPI_Comm_rank")
  call MPI_Comm_size(MPI_COMM_WORLD, info%mpi%nproc, err); call MPI_check(err, "MPI_Comm_size")
#else
  info%mpi%enabled = .false.
  info%mpi%rank    = 0
  info%mpi%nproc   = 1
#endif

! suppress printing for higher ranks
if(info%mpi%rank > 0) isPrint=.false.

! ----- parse command line arguments ----------------------------------------------------

call parse_command_args(cli_opts,err,message)
if(err/=0) stop trim(message)

if(isPrint)then
  print*, 'Control file = ', cli_opts%control_file
  print*, 'Run mode     = ', cli_opts%runmode
endif

! ----- initialize the model domain -----------------------------------------------------

! read hydromet metadata (space/time/coords), apply MPI decomposition, and allocate domain arrays
call setup_domain(cli_opts, info, domain, err, message) 
if(err/=0) stop trim(message)

! ----- initialize model configurations -------------------------------------------------

! choose model, load parameter metadata, derive parameters, and define NetCDF output files
call setup_model_definition(cli_opts, info, work, domain, APAR, BL, BU, err, message)
if(err/=0) stop trim(message)

! ----- initialize work structures ------------------------------------------------------

! allocate space for work structures that depend on number of states and parameters
call init_fuse_work(info, work, err, message)
if(err/=0) stop trim(message)

! ----- set initial counters ------------------------------------------------------------

! Define output and parameter files
ONEMOD=1                      ! one file per model (i.e., model dimension = 1)
work%run%n_evaluations = 0    ! counter for parameter sets evaluated

! ---------------------------------------------------------------------------------------
! ----- run different FUSE modes --------------------------------------------------------
! ---------------------------------------------------------------------------------------

! select fuse mode
select case(trim(cli_opts%runmode))

  ! ----- single parameter set ----------------------------------------------------------

  case('def', 'idx', 'opt')

    OUTPUT_FLAG=.TRUE.

    ! load specific parameter set given index in vector into APAR
    if (cli_opts%runmode=='idx') then
     CALL GET_PRE_PARAM(cli_opts%sets_file, cli_opts%indx, ONEMOD, info%config%nParam, APAR)
    endif

    ! load best parameter set from NetCDF file into APAR
    if (cli_opts%runmode=='opt') then
     CALL GET_SCE_PARAM(cli_opts%sets_file, ONEMOD, info%config%nParam, APAR)
    endif

    ! run FUSE
    CALL FUSE_EVALUATE(APAR, info, work, domain, OUTPUT_FLAG, METRIC_VAL, err, message)
    if(err/=0) stop trim(message)

  ! ----- SCE calibration run -----------------------------------------------------------

  case('sce')

    call sce_driver(info, work, domain, APAR, BL, BU)

  case default
    stop "cannot identify FUSE mode"
  
end select ! (FUSE mode) 
  
! ----- finalize ------------------------------------------------------------------------
  
! deallocate space
DEALLOCATE(APAR, BL, BU, stat=err)
if(err/=0)then; write(*,*) 'unable to deallocate space for parameter vectors'; stop; endif

! close NetCDF files
PRINT *, 'Closing hydromet file'
err = nf90_close(info%files%ncid_hydromet)
if(err/=0)then; message=trim(message)//' nf90_close failed: '//trim(nf90_strerror(err)); return; endif

PRINT *, 'Closing output file'
err = nf90_close(ncid_out)
if(err/=0)then; message=trim(message)//' nf90_close failed: '//trim(nf90_strerror(err)); return; endif

PRINT *, 'Done'
STOP

! ---------------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------

contains

! ----- MPI checker ---------------------------------------------------------------------

subroutine mpi_check(ierr, callee)
#ifdef __MPI__
  use mpi
#endif
  implicit none
  integer(i4b), intent(in) :: ierr
  character(len=*), intent(in) :: callee
#ifdef __MPI__
  integer(i4b) :: slen, ierr2
  character(len=256) :: errstr
  if (ierr /= MPI_SUCCESS) then
    call MPI_Error_string(ierr, errstr, slen, ierr2)
    write(*,*) "MPI error at ", trim(callee), ": ", trim(errstr(1:slen))
    call MPI_Abort(MPI_COMM_WORLD, ierr, ierr2)
  end if
#else
  ! serial build: do nothing
#endif
end subroutine mpi_check

END PROGRAM DISTRIBUTED_DRIVER
