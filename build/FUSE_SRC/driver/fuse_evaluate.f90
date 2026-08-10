MODULE fuse_evaluate_module

  use nrtype
  use multi_flux_types, only: fluxes
  use info_types, only: fuse_info
  use work_types, only: fuse_work
  use domain_types, only: domain_data

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
    ! Modified by Cyril Thébault to include interception, 7/2026
    ! ---------------------------------------------------------------------------------------
    ! Purpose:
    ! --------
    ! Calculate the metric chosen as objective function for single FUSE model and single parameter set
    !   input: model parameter set
    !   output: metric chosen as objective function
    ! ---------------------------------------------------------------------------------------

    use nrtype
    use fuse_globaldata,only: NPAR_SNOW, isPrint, nFUSE_eval
    use multiforce,only: nspat1, nspat2, numtim_sub
    use multibands,only: N_BANDS, n_bands

    IMPLICIT NONE

    ! input
    REAL(WP),DIMENSION(:) , intent(in)     :: XPAR           ! model parameter set
    type(fuse_info)       , intent(in)     :: info           ! info structures (runtime settings etc.)
    type(fuse_work)       , intent(inout)  :: work           ! work structures that depend on npar/nState
    type(domain_data)     , intent(inout)  :: domain         ! the fuse domain structure that stores data arrays
    LOGICAL(LGT)          , intent(in)     :: OUTPUT_FLAG    ! .TRUE. if desire time series output

    ! output
    REAL(WP),INTENT(OUT)                   :: METRIC_VAL     ! metric 

    ! error control
    integer(i4b)                           :: err, ierr
    character(len=1024)                    :: message

    ! timing
    real(wp)                               :: t1, t2

    ! ---------------------------------------------------------------------------------------

    ! populate parameter structures and initialize states
    call initialize_run(info, XPAR, work, domain, ierr, message)
    if (ierr /= 0) stop trim(message)
    
    ! initialize timing
    CALL CPU_TIME(T1)

    ! run fuse for the entire time series
    call run_time_loop(info, work, domain, OUTPUT_FLAG, ierr, message)
    if (ierr /= 0) stop trim(message)

    ! get timing information
    CALL CPU_TIME(T2)
    if(isPrint) WRITE(*,*) "TIME ELAPSED = ", t2-t1

    ! calculate mean summary statistics
    if (info%mpi%rank == 0) then

      if(isPrint) PRINT *, 'Calculating performance metrics...'
      CALL MEAN_STATS(work, domain)
      metric_val = work%run%stats%metric_val

      write(*,'(i6,1x,a11,1x,f12.6,1x,a20,1x,f12.6)') nFUSE_eval, "OBJ FUNC = ", METRIC_VAL, "; TIME ELAPSED = ", t2-t1
      !if(nFUSE_eval > 10) stop "checking results"

      if(isPrint) PRINT *, 'Writing model statistics...'
      CALL PUT_SSTATS(work%run%stats, work%run%n_evaluations)

    endif ! all observations on rank=0

  END SUBROUTINE fuse_evaluate

  ! -------------------------------------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------------------------------------

  ! -------------------------------------------------------------------------------------------------------------------
  ! ----- private subroutine initialize_run: populate param sets and initialize states  -------------------------------
  ! -------------------------------------------------------------------------------------------------------------------

  subroutine initialize_run(info, xpar, work, domain, ierr, message)
  
  use fuse_globaldata,  only: isPrint, fracstate0
  use model_defn,  only: SMODL
  use model_defnames
  
  use multiforce,  only: nspat1, nspat2
  use multistate,  only: FSTATE
  use multibands
 
  use par_insert_module, only: put_parset
  use par_derive_module, only: par_derive 
  use par_insert_module
  use str_2_xtry_module
  use xtry_2_str_module
  use put_params_module, only: put_params

  use qtimedelay_module, only: qtimedelay
  use init_state_module, only: init_state

  implicit none

  type(fuse_info)        , intent(in)      :: info
  real(wp), dimension(:) , intent(in)      :: xpar
  type(fuse_work)        , intent(inout)   :: work

  type(domain_data)      , intent(inout)   :: domain
  integer(i4b)           , intent(out)     :: ierr
  character(len=*)       , intent(out)     :: message

  integer(i4b)                             :: iSpat1, iSpat2, iBands
  character(len=256)                       :: cmessage  ! error message of downwind routine

  ierr    = 0
  message = "initialize_run/"

  ! increment parameter counter for model output
  work%run%n_evaluations = work%run%n_evaluations + 1

  ! add parameter set to the data structure
  call put_parset(xpar, info%config%listParam, work%par, ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  if (isPrint) then
    print *, 'Parameter set added to data structure:'
    print *, xpar
  end if

  ! Compute derived model parameters (bucket sizes, etc.)
  CALL PAR_DERIVE(work%par)   ! populates work structure with DPARAM

  ! compute fraction of runoff in future time steps (unit hydrograph for hillslope routing)
  CALL QTIMEDELAY(info,                  & !
                  work%par%param_adjust, & ! adjustable model parameters (time delay, )
                  work%par%param_derive, & ! derived model parameters (FRAC_FUTURE, )
                  work%route%future,     & ! rolling routing convolution queue
                  ierr, cmessage)
  if (ierr/=0) then; message=trim(message)//trim(cmessage); return; endif

  ! get elevation bands (if catchment)
  Z_FORCING      = domain%z_forcing(1,1)
  MBANDS(:)%info = domain%bands_info(1,1,:)

  if (isPrint) print *, 'Writing parameter values...'
  call put_params(work%run%n_evaluations)

  ! initialize model states over the 2D gridded domain (1 x nHRU in catchment mode)
  do iSpat2 = 1, nSpat2
    do iSpat1 = 1, nSpat1
      call init_state(fracstate0,            & ! input:  fraction state
                      work%par%param_adjust, & ! adjustable model parameters (time delay, )
                      work%par%param_derive, & ! derived model parameters (FRAC_FUTURE, )
                      work%step%state0,      & ! output: start-of-step state
                      work%snow%sbands(:)%var%bands_var)  ! output: SWE for elevation bands
      call str_2_xtry(work%step%state0, work%num%x0)
      domain%state(iSpat1, iSpat2, 1) = work%step%state0
    end do
  end do
  if (isPrint) print *, 'Model states initialized over the 2D gridded domain'

  ! initialize elevation bands if snow module is on
  if (isPrint) print *, 'N_BANDS =', N_BANDS
  if (SMODL%iSNOWM == iopt_temp_index) then

    ! initialize template once (spat1, spat2, bands, first time index)
    domain%bands_var(:,:,:,1)%SWE         = 0._wp
    domain%bands_var(:,:,:,1)%SNOWACCMLTN = 0._wp
    domain%bands_var(:,:,:,1)%SNOWMELT    = 0._wp
    domain%bands_var(:,:,:,1)%DSWE_DT     = 0._wp
   
    ! get work array foor the bands
    work%snow%sbands(:)%var%bands_var = domain%bands_var(1,1,:,1)

  end if ! if (SMODL%iSNOWM == iopt_temp_index)

  ! initialize summary statistics + timer
  call init_stats(work)

  if (isPrint) print *, 'End of initialize_run'

  end subroutine initialize_run

  ! -------------------------------------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------------------------------------

  ! -------------------------------------------------------------------------------------------------------------------
  ! ----- private subroutine run_time_loop: run fuse for the entire time series  --------------------------------------
  ! -------------------------------------------------------------------------------------------------------------------

  subroutine run_time_loop(info, work, domain, output_flag, ierr, message)

  use fuse_globaldata, only: isPrint
  use multiforce, only: timDat  ! NOTE: used in legacy codes
  use multiforce, only: nspat1, nspat2, DELTIM, sim_beg, sim_end, numtim_sub
  use time_utils,        only: caldatss
  use get_hydromet_module, only: get_met_data
  use get_hydromet_module, only: get_qobs_data
  use put_output_module,   only: put_output

  implicit none

  type(fuse_info)   , intent(in)    :: info           ! info structures that include "everything"
  type(fuse_work)   , intent(inout) :: work           ! work structures that depend on npar/nState
  type(domain_data) , intent(inout) :: domain         ! domain structures that hold 3-d data 
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
  real(wp)     :: dt_sub, dt_full
  integer(i4b) :: iSpat1, iSpat2, iBands
  character(len=256) :: cmessage  

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
    ! ----- start of subperiod: load hydromet data --------------------------------------------------------------------

    ! determine length of current subperiod
    remaining = sim_end - in_idx + 1       ! # remaining data windows in simulation 
    chunk_len = min(numtim_sub, remaining) ! # data windows in the sub-period
 
    ! save the start of the chunks (avoid arithmetic)
    chunk_start_sim = in_idx - sim_beg + 1     ! start of chunk in simulation index space
    chunk_start_in  = in_idx                   ! start of chunk in input index space

    if(isPrint) PRINT *, 'New subperiod: loading hydromet data for ',chunk_len,' time steps'
    
    ! load meteorological forcing data for desired period into the domain%force data structure
    ! NOTE: reads different spatial slice per MPI rank in the nSpat2 dimension
    call get_met_data(info, chunk_start_in, chunk_len, &
                      domain, ierr, cmessage)
    if (ierr/=0) then; message=trim(message)//trim(cmessage); return; endif

    ! load streamflow observations for desired period into the domain%valid data structure
    ! NOTE: qobs replicated across MPI ranks
    call get_qobs_data(info, chunk_start_in, chunk_len, &
                       work%obs%q(:,1:chunk_len), ierr, cmessage)
    if (ierr/=0) then; message=trim(message)//trim(cmessage); return; endif

    if(isPrint) PRINT *, 'Hydromet data loaded. Running FUSE...'
    
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

      ! loop through grid points and run the model for one time step
      DO iSpat2=1,nSpat2
        DO iSpat1=1,nSpat1
  
          ! run fuse for one grid cell
          call advance_one_cell(work, domain, sub_idx, iSpat1, iSpat2, dt_sub, dt_full, ierr, message)
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
      CALL PUT_OUTPUT(info, work, domain, chunk_start_sim, chunk_start_in, chunk_len)
      if(isPrint) PRINT *, 'Done writing output'
    ELSE
      if(isPrint) PRINT *, 'OUTPUT_FLAG is set on FALSE, no output written'
    END IF

    ! TODO: set domain%state and domain%bands_var to NA

    ! reinitialize states for next subperiod using last time step
    domain%state(:,:,1)       = domain%state(:,:,chunk_len+1)
    domain%bands_var(:,:,:,1) = domain%bands_var(:,:,:,chunk_len+1)

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

  subroutine advance_one_cell(work, domain, sub_idx, iSpat1, iSpat2, dt_sub, dt_full, ierr, message)

  ! switches / options
  use fuse_globaldata,   only: NA_VALUE_SP
  use set_all_module, only: SET_STATE, SET_FLUXES, SET_ROUTE

  ! model options
  use model_defn,   only: SMODL
  use model_defnames

  ! physics drivers
  use physics_orig_module, only: physics_orig
  use physics_diff_module, only: physics_diff

  use Q_OVERLAND_module, only: Q_OVERLAND

  ! state vector conversions
  use str_2_xtry_module
  use xtry_2_str_module

  ! diff-mode flags
  use model_numerix, only: diff_mode, original, differentiable

  implicit none

  type(fuse_work)       , intent(inout) :: work           ! work structures that depend on npar/nState
  type(domain_data)     , intent(inout) :: domain         ! domain structures that hold 3-d data 
  integer(i4b)          , intent(in)    :: sub_idx, iSpat1, iSpat2
  real(wp)              , intent(inout) :: dt_sub, dt_full
  integer(i4b)          , intent(out)   :: ierr
  character(len=*)      , intent(out)   :: message

  ! locals
  character(len=1024)       :: cmessage

  ierr = 0
  message = "advance_one_cell/"

  ! only run FUSE for grid points within domain when elev_mask is FALSE
  if (.not. domain%elev_mask(iSpat1,iSpat2)) then

    ! -----------------------------------------------------------------------------------
    ! ----- initialize ------------------------------------------------------------------
    ! -----------------------------------------------------------------------------------

    ! extract forcing for this grid cell and time step
    work%step%force = domain%force(iSpat1,iSpat2,sub_idx)

    call check_force(work%step%force%ppt, work%step%force%temp, ierr, cmessage)
    if (ierr /= 0) then; message = trim(message)//trim(cmessage); return; end if

    ! extract model states for this grid cell and time step
    work%step%state0 = domain%state(iSpat1,iSpat2,sub_idx)
    work%step%state1 = domain%state(iSpat1,iSpat2,sub_idx)
    call STR_2_XTRY(work%step%state0, work%num%x0)

    ! extract the elevation bands
    work%snow%z_forcing               = domain%z_forcing (iSpat1,iSpat2)
    work%snow%sbands(:)%info          = domain%bands_info(iSpat1,iSpat2,:)
    work%snow%sbands(:)%var%bands_var = domain%bands_var (iSpat1,iSpat2,:,sub_idx)

    ! initialize model fluxes
    call INITFLUXES(work%step%flux)

    ! initialize the snow fluxes for the elevation bands
    work%snow%sbands(:)%var%SNOWACCMLTN = 0._wp
    work%snow%sbands(:)%var%SNOWMELT    = 0._wp

    ! initialize derivatives in the adjoint model
    work%adj%df_dS(:)   = work%step%flux
    work%adj%df_dPar(:) = work%step%flux

    ! -----------------------------------------------------------------------------------
    ! ----- run -------------------------------------------------------------------------
    ! -----------------------------------------------------------------------------------

    select case (diff_mode)

      case (original)
        call physics_orig(work,                     &
                          dt_sub, dt_full,          &
                          sub_idx, iSpat1, iSpat2,  &
                          ierr, cmessage)

      case (differentiable)
        call physics_diff(work, dt_full,            &
                          sub_idx, iSpat1, iSpat2,  &
                          ierr, cmessage)

      case default
        ierr     = 10
        cmessage = 'unknown physics mode'

    end select

    if (ierr /= 0) then
      message = trim(message)//trim(cmessage)
      return
    end if

    ! -----------------------------------------------------------------------------------
    ! ----- routing ---------------------------------------------------------------------
    ! -----------------------------------------------------------------------------------
   
    ! hillslope routing
    call Q_OVERLAND(work%par%param_derive, & ! params
                    work%step%flux,        & ! land fluxes
                    work%route%future,     & ! rolling routing convolution queue
                    work%step%route,       & ! routing fluxes
                    ierr, cmessage)
    if (ierr /= 0) then; message = trim(message)//trim(cmessage); return; end if

    call check_routing(work%step%route%q_routed, ierr, cmessage)
    if (ierr /= 0) then; message = trim(message)//trim(cmessage); return; end if

    ! -----------------------------------------------------------------------------------
    ! ----- finalize ---------------------------------------------------------------------
    ! -----------------------------------------------------------------------------------
    
    ! convert to data structures
    call XTRY_2_STR(work%num%x1, work%step%state1)

    ! write back to 3D data structures
    domain%state(iSpat1,iSpat2,sub_idx+1) = work%step%state1
    domain%flux (iSpat1,iSpat2,sub_idx)   = work%step%flux
    domain%route(iSpat1,iSpat2,sub_idx)   = work%step%route

    ! save canonical snow state
    if (SMODL%iSNOWM == iopt_temp_index) then

      domain%state(iSpat1,iSpat2,sub_idx+1)%SWE_TOT = sum(work%snow%sbands(:)%var%SWE * &
                                                          work%snow%sbands(:)%info%AF)

      domain%bands_var(iSpat1,iSpat2,:,sub_idx+1)   = work%snow%sbands(:)%var%bands_var

    end if

    ! stats
    call COMP_STATS(work)

  ! -------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------
  
  ! outside mask: NA fill
  
  else
    
    call SET_STATE (NA_VALUE_SP, domain%state(iSpat1,iSpat2,sub_idx+1), &
                                 domain%bands_var(iSpat1,iSpat2,:,sub_idx+1) )

    call SET_FLUXES(NA_VALUE_SP, domain%flux(iSpat1,iSpat2,sub_idx),    &
                                 domain%bands_var(iSpat1,iSpat2,:,sub_idx+1) )

    call SET_ROUTE (NA_VALUE_SP, domain%route(iSpat1,iSpat2,sub_idx) )
  
  end if

  ! -------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------

  contains

    ! -----------------------------------------------------------------------------------
    ! -----------------------------------------------------------------------------------
    
    subroutine check_force(ppt, temp, ierr, message)
    
    real(wp)               , intent(in)    :: ppt
    real(wp)               , intent(in)    :: temp
      
    integer(i4b)           , intent(out)   :: ierr
    character(*)           , intent(out)   :: message

    ierr    = 0
    message = "check_force/"

    ! forcing sanity checks
    if (ppt < 0.0_wp) then
      ierr=1; message='Negative precipitation in input file'; return
    end if

    if (ppt > 5000.0_wp) then
      ierr=1; message='Precipitation greater than 5000 in input file'; return
    end if

    if (ppt < 0.0_wp) then
      ierr=1; message='Negative PET in input file'; return
    end if

    if (ppt > 100.0_wp) then
      ierr=1; message='PET greater than 100 in input file'; return
    end if

    if (temp < -100.0_wp) then
      ierr=1; message='Temperature lower than -100 in input file'; return
    end if

    if (temp > 100.0_wp) then
      ierr=1; message='Temperature greater than 100 in input file'; return
    end if

    end  subroutine check_force

    ! -----------------------------------------------------------------------------------
    ! -----------------------------------------------------------------------------------
   
    subroutine check_routing(q_routed, ierr, message) 
    
    real(wp)               , intent(in)    :: q_routed
      
    integer(i4b)           , intent(out)   :: ierr
    character(*)           , intent(out)   :: message

    ierr    = 0
    message = "check_force/"
    ierr    = 0
    message = "check_routing/"

    ! routing sanity checks

    if (q_routed < 0._wp) then
      message=trim(message)//'Q_ROUTED is less than zero'
      ierr=1; return
    end if

    if (q_routed > 1000._wp) then
      message=trim(message)//'Q_ROUTED is enormous'
      ierr=1; return
    end if

    end subroutine check_routing

    ! -----------------------------------------------------------------------------------
    ! -----------------------------------------------------------------------------------

  end subroutine advance_one_cell

END MODULE fuse_evaluate_module
