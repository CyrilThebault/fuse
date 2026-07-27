MODULE fuse_evaluate_module

  use nrtype
  use multi_flux_types, only: fluxes
  use info_types, only: fuse_info
  use work_types, only: fuse_work
  use data_types, only: domain_data

  IMPLICIT NONE

  CONTAINS
  
    SUBROUTINE fuse_evaluate(XPAR, info, work, domain, OUTPUT_FLAG, METRIC_VAL)

    ! ---------------------------------------------------------------------------------------
    ! Creator:
    ! --------
    ! Martyn Clark, 2009
    ! Modified by Brian Henn to include snow model, 6/2013
    ! Modified by Nans Addor to enable grid-based modeling, 9/2016
    ! Modified by Cyril Thébault to allow different metrics as objective function, 2024
    ! Modified by Martyn Clark to call differentiable modeling routines, 12/2025
    ! Modified by Martyn Clark to simplify/refactor, 02/2026
    ! ---------------------------------------------------------------------------------------
    ! Purpose:
    ! --------
    ! Calculate the metric chosen as objective function for single FUSE model and single parameter set
    !   input: model parameter set
    !   output: metric chosen as objective function
    ! ---------------------------------------------------------------------------------------

    use nrtype
    use fuse_globaldata,only: NPAR_SNOW, isPrint, nFUSE_eval
    use model_defn,only: NSTATE
    use multiparam,only: NUMPAR
    use multiforce,only: nspat1, nspat2, numtim_sub
    use multibands,only: N_BANDS, n_bands
    use multistats,only: MSTATS, PCOUNT
    use multi_flux,only: W_FLUX_3d

    IMPLICIT NONE

    ! input
    REAL(SP),DIMENSION(:) , intent(in)     :: XPAR           ! model parameter set
    type(fuse_info)       , intent(in)     :: info           ! info structures (runtime settings etc.)
    type(fuse_work)       , intent(inout)  :: work           ! work structures that depend on npar/nState
    type(domain_data)     , intent(inout)  :: domain         ! the fuse domain structure that stores data arrays
    LOGICAL(LGT)          , intent(in)     :: OUTPUT_FLAG    ! .TRUE. if desire time series output

    ! output
    REAL(SP),INTENT(OUT)                   :: METRIC_VAL     ! metric 

    ! error control
    integer(i4b)                           :: err, ierr
    character(len=1024)                    :: message

    ! timing
    real(sp)                               :: t1, t2

    ! ---------------------------------------------------------------------------------------

    ! allocate 3d data structure for fluxes
    allocate(w_flux_3d(nspat1, nspat2, numtim_sub), stat=ierr)
    if (ierr /= 0) stop "problem allocating w_flux_3d in fuse_evaluate"

    ! populate parameter structures and initialize states
    call initialize_run(XPAR, work, ierr, message)
    if (ierr /= 0) stop trim(message)
    
    ! initialize timing
    CALL CPU_TIME(T1)

    ! run fuse for the entire time series
    call run_time_loop(info, work, OUTPUT_FLAG, err, message)
    if (err /= 0) stop trim(message)

    ! get timing information
    CALL CPU_TIME(T2)
    if(isPrint) WRITE(*,*) "TIME ELAPSED = ", t2-t1

    ! calculate mean summary statistics
    ! NOTE: .NOT.GRID_FLAG means catchment mode (lumped or distributed)
    if( .not. info%space%grid_flag)then

      if(isPrint) PRINT *, 'Calculating performance metrics...'
      CALL MEAN_STATS()
      METRIC_VAL = MSTATS%METRIC_VAL

      write(*,'(i6,1x,a11,1x,f12.6,1x,a20,1x,f12.6)') nFUSE_eval, "OBJ FUNC = ", METRIC_VAL, "; TIME ELAPSED = ", t2-t1
      !if(nFUSE_eval > 10) stop "checking results"

    endif ! if catchment mode (lumped or distributed)

    if(isPrint) PRINT *, 'Writing model statistics...'
    CALL PUT_SSTATS(PCOUNT)

    ! deallocate output buffer
    DEALLOCATE(W_FLUX_3d); IF (IERR.NE.0) STOP ' problem deallocating W_FLUX_3d in fuse_metric '

  END SUBROUTINE fuse_evaluate

  ! -------------------------------------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------------------------------------

  ! -------------------------------------------------------------------------------------------------------------------
  ! ----- private subroutine initialize_run: populate param sets and initialize states  -------------------------------
  ! -------------------------------------------------------------------------------------------------------------------

  subroutine initialize_run(xpar, work, err, message)
  
  use fuse_globaldata,  only: isPrint, fracstate0
  use model_defn,  only: SMODL
  use model_defnames
  
  use multiparam,  only: NUMPAR
  use multiforce,  only: nspat1, nspat2, DELTIM
  use multistate,  only: FSTATE, gState_3d
  use multistats,  only: PCOUNT
  use multibands
 
  use par_derive_module, only: par_derive 
  use par_insert_module
  use str_2_xtry_module
  use xtry_2_str_module
  use put_params_module, only: put_params
  implicit none

  real(sp), dimension(:) , intent(in)      :: xpar
  type(fuse_work)        , intent(inout)   :: work

  integer(i4b)           , intent(out)     :: err
  character(len=*)       , intent(out)     :: message

  integer(i4b) :: iSpat1, iSpat2, iBands

  err = 0
  message = ""

  ! increment parameter counter for model output
  PCOUNT = PCOUNT + 1

  ! add parameter set to the data structure
  call put_parset(xpar)
  if (isPrint) then
    print *, 'Parameter set added to data structure:'
    print *, xpar
  end if

  ! compute derived model parameters (bucket sizes, etc.)
  call par_derive(err, message)
  if (err /= 0) then
    write(*,*) trim(message)
    stop
  end if

  ! get elevation bands (if catchment)
  Z_FORCING      = Z_FORCING_grid(1,1)
  MBANDS(:)%info = MBANDS_INFO_3d(1,1,:)

  if (isPrint) print *, 'Writing parameter values...'
  call put_params(PCOUNT)

  ! initialize model states over the 2D gridded domain (1x1 in catchment mode)
  do iSpat2 = 1, nSpat2
    do iSpat1 = 1, nSpat1
      call init_state(fracstate0)
      call str_2_xtry(FSTATE, work%num%x0)
      call xtry_2_str(work%num%x0, FSTATE)
      gState_3d(iSpat1, iSpat2, 1) = FSTATE
    end do
  end do
  if (isPrint) print *, 'Model states initialized over the 2D gridded domain'

  ! initialize elevation bands if snow module is on
  if (isPrint) print *, 'N_BANDS =', N_BANDS
  if (SMODL%iSNOWM == iopt_temp_index) then

    ! initialize template once
    work%snow%sbands(:)%var%SWE         = 0._sp
    work%snow%sbands(:)%var%SNOWACCMLTN = 0._sp
    work%snow%sbands(:)%var%SNOWMELT    = 0._sp
    work%snow%sbands(:)%var%DSWE_DT     = 0._sp

    ! copy to every grid cell (legacy staging)
    do iSpat2 = 1, nSpat2
      do iSpat1 = 1, nSpat1
        do iBands = 1, n_bands
          MBANDS_VAR_4d(iSpat1, iSpat2, iBands, 1) = work%snow%sbands(iBands)%var%bands_var
        end do
      end do
    end do

    if (isPrint) print *, 'Snow states initialized over the 2D gridded domain'
  end if

  ! initialize summary statistics + timer
  call init_stats()

  end subroutine initialize_run

  ! -------------------------------------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------------------------------------

  ! -------------------------------------------------------------------------------------------------------------------
  ! ----- private subroutine run_time_loop: run fuse for the entire time series  --------------------------------------
  ! -------------------------------------------------------------------------------------------------------------------

  subroutine run_time_loop(info, work, output_flag, ierr, message)

  use fuse_globaldata, only: isPrint
  use multiforce, only: timDat  ! NOTE: used in legacy cides
  use multiforce, only: nspat1, nspat2, DELTIM, sim_beg, sim_end, numtim_sub
  use multistate, only: gState_3d
  use multibands, only: MBANDS_VAR_4d
  use time_utils,        only: caldatss
  use getPETgrid_module, only: getPETgrid
  use get_gforce_module, only: get_gforce_3d
  use put_output_module, only: put_output

  implicit none

  type(fuse_info)   , intent(in)    :: info           ! info structures that include "everything"
  type(fuse_work)   , intent(inout) :: work           ! work structures that depend on npar/nState
  logical(lgt)      , intent(in)    :: output_flag

  integer(i4b)      , intent(out)   :: ierr
  character(len=*)  , intent(out)   :: message

  ! time management
  integer(i4b) :: sim_idx           ! index of simulation: 1..numtim_sim
  integer(i4b) :: sub_idx           ! index of forcing slice: 1..chunk_len
  integer(i4b) :: in_idx            ! index of input NetCDF time axis: sim_beg..sim_end
  integer(i4b) :: remaining         ! # remaining data windows in simulation
  integer(i4b) :: chunk_len         ! # data windows in the sub-period
  integer(i4b) :: chunk_start_in    ! start-of-chunk index in the input file
  integer(i4b) :: chunk_start_sim   ! start-of-chunk index in the simulation

  ! locals
  logical(lgt), parameter :: computePET = .false.
  real(sp)     :: dt_sub, dt_full
  integer(i4b) :: iSpat1, iSpat2, iBands

  ierr = 0
  message = "run_time_loop/"

  ! This version of FUSE enables the user to load slices of the forcing
  !
  ! FUSE1 used to access the input file at each time step, slowing operations
  ! down over large domains on systems with slow I/O. The number of timesteps
  ! of the slices is defined by the user in the filemanager. The default is
  ! that the whole time period needed for the simulation is loaded, but
  ! this can exceed memory capacity when large domains are processed.
  
  ! To overcome this, a subperiod (slice) of the forcing can be loaded in
  ! memory and used to run FUSE. Then, the results are saved to the
  ! output file, and the next slice of forcing is loaded. This enables FUSE to
  ! run quicker than when forcing is loaded at each time step and grid point,
  ! while also controlling memory usage.

  ! initialize model time step
  dt_sub  = DELTIM
  dt_full = DELTIM

  ! initialise time indices for whole simulation and subperiod
  sub_idx = 1 ! index in data subset

  ! ----- loop through chunks of data ---------------------------------------------------------------------------------

  in_idx = sim_beg
  do while (in_idx <= sim_end)

    ! get the simulation index
    sim_idx = in_idx - sim_beg + 1

    ! -----------------------------------------------------------------------------------------------------------------
    ! ----- start of subperiod: load forcing --------------------------------------------------------------------------

    ! determine length of current subperiod
    remaining = sim_end - in_idx + 1       ! # remaining data windows in simulation 
    chunk_len = min(numtim_sub, remaining) ! # data windows in the sub-period
 
    ! save the start of the chunks (avoid arithmetic)
    chunk_start_sim = in_idx - sim_beg + 1     ! start of chunk in simulation index space
    chunk_start_in  = in_idx                   ! start of chunk in input index space

    ! load forcing for desired period into gForce_3d
    if(isPrint) PRINT *, 'New subperiod: loading forcing for ',chunk_len,' time steps'
    call get_gforce_3d(info, chunk_start_in, chunk_len, ierr, message)
    IF(ierr/=0) stop 'Error while extracting 3d forcing: '//trim(message)
    if(isPrint) PRINT *, 'Forcing loaded. Running FUSE...'
    
    ! -----------------------------------------------------------------------------------------------------------------

    ! -----------------------------------------------------------------------------------------------------------------
    ! ----- loop through data chunk (sub-period) ----------------------------------------------------------------------

    do sub_idx = 1, chunk_len

      ! get indices in the input file (in_idx) and the simulation period (sim_idx)
      in_idx  = chunk_start_in + sub_idx - 1
      sim_idx = chunk_start_sim + sub_idx - 1

      ! get the model time
      call caldatss(info%time%jdate(in_idx), work%step%time%iy, work%step%time%im, work%step%time%id, &
                                             work%step%time%ih, work%step%time%imin, work%step%time%dsec)
      timDat = work%step%time ! NOTE: used in the legacy data structures

      ! compute potential ET
      IF(computePET) CALL getPETgrid(ierr,message)
      IF(ierr/=0) stop TRIM(message)
   
      ! loop through grid points and run the model for one time step
      DO iSpat2=1,nSpat2
        DO iSpat1=1,nSpat1
  
          ! run fuse for one grid cell
          call advance_one_cell(work, sub_idx, iSpat1, iSpat2, dt_sub, dt_full, ierr, message)
          if (ierr /= 0)  stop trim(message)
  
          !if(sub_idx > 100) stop "check"

        END DO  ! (looping thru 2nd spatial dimension)
      END DO  ! (looping thru 1st spatial dimension)

    end do  ! looping through subperiod

    !stop "looping through time period"

    ! -----------------------------------------------------------------------------------------------------------------
    
    ! -----------------------------------------------------------------------------------------------------------------
    ! ----- end of subperiod: write to output file and save states ----------------------------------------------------

    if(isPrint) PRINT *, 'End of subperiod reached:'

    ! write model output
    IF (OUTPUT_FLAG) THEN
      if(isPrint) PRINT *, 'Write output for ',chunk_len,' time steps starting at indices', chunk_start_sim
      CALL PUT_OUTPUT(work, chunk_start_sim, chunk_start_in, chunk_len)
      if(isPrint) PRINT *, 'Done writing output'
    ELSE
      if(isPrint) PRINT *, 'OUTPUT_FLAG is set on FALSE, no output written'
    END IF

    ! TODO: set gState_3d and MBANDS_VAR_4d to NA

    ! reinitialize states for next subperiod using last time step
    gState_3d(:,:,1)       = gState_3d(:,:,chunk_len+1)
    MBANDS_VAR_4d(:,:,:,1) = MBANDS_VAR_4d(:,:,:,chunk_len+1)

    ! -----------------------------------------------------------------------------------------------------------------

    ! update the index in the input file
    in_idx = chunk_start_in + chunk_len

  END DO  ! (loop through timesteps)

  end subroutine run_time_loop
  
  ! -------------------------------------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------------------------------------

  ! -------------------------------------------------------------------------------------------------------------------
  ! ----- private subroutine advance_one_cell: run fuse for one grid cell ---------------------------------------------
  ! -------------------------------------------------------------------------------------------------------------------

  subroutine advance_one_cell(work, sub_idx, iSpat1, iSpat2, dt_sub, dt_full, err, message)

  ! switches / options
  use fuse_globaldata,   only: NA_VALUE_SP
  use model_defn,   only: SMODL, NSTATE
  use model_defnames
  use multiforce,   only: DELTIM, gForce_3d, aForce, MFORCE, nspat1, nspat2
  use multistate,   only: gState_3d, FSTATE, MSTATE
  use multiroute,   only: MROUTE, AROUTE_3d
  use multibands
  use multi_flux,   only: W_FLUX, W_FLUX_3d
  use set_all_module, only: SET_STATE, SET_FLUXES, SET_ROUTE

  ! state vector conversions
  use str_2_xtry_module
  use xtry_2_str_module

  ! differentiable
  use get_bundle_module,       only: get_bundle
  use implicit_solve_module,   only: implicit_solve
  use update_swe_diff_module,  only: update_swe_diff  ! (only if you actually call it here)
  use update_swe_diff_module,  only: update_swe_diff  ! ok to remove if unused

  ! original solver interface
  use interfaceb, only: ode_int, fuse_solve

  ! diff-mode flags (make sure these names really live here in your tree)
  use model_numerix, only: diff_mode, original, differentiable

  implicit none

  type(fuse_work)       , intent(inout) :: work           ! work structures that depend on npar/nState
  integer(i4b)          , intent(in)    :: sub_idx, iSpat1, iSpat2
  real(sp)              , intent(inout) :: dt_sub, dt_full
  integer(i4b)          , intent(out)   :: err
  character(len=*)      , intent(out)   :: message

  ! locals
  integer(i4b)              :: ierr
  character(len=1024)       :: cmessage

  err = 0
  message = "advance_one_cell/"
  ierr = 0
  cmessage = ""

  ! ---------------------------------------------------------------------------
  ! only run FUSE for grid points within domain defined by elev_mask
  ! NOTE: you currently run when elev_mask is FALSE (keep as-is for BFB)
  ! ---------------------------------------------------------------------------
  if (.not. elev_mask(iSpat1,iSpat2)) then

    ! extract forcing for this grid cell and time step
    MFORCE = gForce_3d(iSpat1,iSpat2,sub_idx)

    ! forcing sanity checks (keep behavior; convert STOP -> error return)
    if (MFORCE%PPT < 0.0_sp) then
      err=1; message='Negative precipitation in input file'; return
    end if
    if (MFORCE%PPT > 5000.0_sp) then
      err=1; message='Precipitation greater than 5000 in input file'; return
    end if
    if (MFORCE%PET < 0.0_sp) then
      err=1; message='Negative PET in input file'; return
    end if
    if (MFORCE%PET > 100.0_sp) then
      err=1; message='PET greater than 100 in input file'; return
    end if
    if (MFORCE%TEMP < -100.0_sp) then
      err=1; message='Temperature lower than -100 in input file'; return
    end if
    if (MFORCE%TEMP > 100.0_sp) then
      err=1; message='Temperature greater than 100 in input file'; return
    end if

    ! extract model states for this grid cell and time step
    FSTATE = gState_3d(iSpat1,iSpat2,sub_idx)
    MSTATE = FSTATE
    call STR_2_XTRY(FSTATE, work%num%x0)

    ! initialize model fluxes
    ! If INITFLUXES lives somewhere else in your tree, swap this line accordingly.
    call INITFLUXES()

    ! populate fuse work structure (diff path only)
    if (diff_mode == differentiable) call get_bundle(work)

    ! -------------------------
    ! snow module
    ! -------------------------
    select case(SMODL%iSNOWM)

      case(iopt_temp_index)

        Z_FORCING      = Z_FORCING_grid(iSpat1,iSpat2)
        MBANDS(:)%info = MBANDS_INFO_3d(iSpat1,iSpat2,:)
        MBANDS(:)%var  = MBANDS_VAR_4d(iSpat1,iSpat2,:,sub_idx)

        if (diff_mode == differentiable) then
          work%snow%z_forcing               = Z_FORCING
          work%snow%sbands(:)%info          = MBANDS(:)%info
          work%snow%sbands(:)%var%bands_var = MBANDS(:)%var
        end if

        select case(diff_mode)
          case(original)
            call UPDATE_SWE(DELTIM)
          case(differentiable)
            call UPDATE_SWE_DIFF(work, DELTIM)
          case default
            err=1; message='advance_one_cell: cannot identify diff_mode (snow)'; return
        end select

      case(iopt_no_snowmod)
        call QRAINERROR()

      case default
        err=1; message='advance_one_cell: unknown SMODL%iSNOWM option'; return

    end select

    ! -------------------------
    ! soil physics
    ! -------------------------
    select case(diff_mode)

      case(original)
        call ODE_INT(FUSE_SOLVE, work%num%x0, work%num%x1, dt_sub, dt_full, ierr, cmessage)
        if (ierr /= 0) then
          err=1; message=trim(cmessage); return
        end if
        !print*, 'original = ', mstate%watr_1, mstate%watr_2, w_flux%QSURF, w_flux%QBASE_2

      case(differentiable)
        call implicit_solve(work, work%num%x0, work%num%x1, nState, ierr, cmessage)
        if (ierr /= 0) then
          err=1; message=trim(cmessage); return
        end if
        W_FLUX = work%step%flux

      case default
        err=1; message='advance_one_cell: cannot identify diff_mode (soil)'; return

    end select

    ! routing
    call Q_OVERLAND()
    if (MROUTE%Q_ROUTED < 0._sp) then
      err=1; message='Q_ROUTED is less than zero'; return
    end if
    if (MROUTE%Q_ROUTED > 1000._sp) then
      err=1; message='Q_ROUTED is enormous'; return
    end if

    ! write back to 3D buffers
    call XTRY_2_STR(work%num%x1, FSTATE)
    gState_3d(iSpat1,iSpat2,sub_idx+1) = FSTATE
    W_FLUX_3d(iSpat1,iSpat2,sub_idx)   = W_FLUX
    AROUTE_3d(iSpat1,iSpat2,sub_idx)   = MROUTE

    if (SMODL%iSNOWM == iopt_temp_index) then

      if (diff_mode == differentiable) then
        Z_FORCING      = work%snow%z_forcing
        MBANDS(:)%info = work%snow%sbands(:)%info
        MBANDS(:)%var  = work%snow%sbands(:)%var%bands_var
      end if

      gState_3d(iSpat1,iSpat2,sub_idx+1)%SWE_TOT = sum(MBANDS(:)%var%SWE * MBANDS(:)%info%AF)
      MBANDS_VAR_4d(iSpat1,iSpat2,:,sub_idx+1)   = MBANDS(:)%var

    end if

    ! forcing diagnostics
    aForce(sub_idx)%ppt = sum(gForce_3d(:,:,sub_idx)%ppt) / real(size(gForce_3d(:,:,sub_idx)), kind=sp)
    aForce(sub_idx)%pet = sum(gForce_3d(:,:,sub_idx)%pet) / real(size(gForce_3d(:,:,sub_idx)), kind=sp)

    ! stats
    call COMP_STATS()

  else
    ! outside mask: NA fill
    call SET_STATE(NA_VALUE_SP)
    gState_3d(iSpat1,iSpat2,sub_idx) = FSTATE

    call SET_FLUXES(NA_VALUE_SP)
    W_FLUX_3d(iSpat1,iSpat2,sub_idx) = W_FLUX

    call SET_ROUTE(NA_VALUE_SP)
    AROUTE_3d(iSpat1,iSpat2,sub_idx) = MROUTE
  end if

  end subroutine advance_one_cell

END MODULE fuse_evaluate_module
