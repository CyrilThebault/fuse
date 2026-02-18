MODULE fuse_evaluate_module

  use nrtype
  use multi_flux_types, only: fluxes
  use work_types, only: fuse_work

  IMPLICIT NONE

  ! temporary type: run context
  type :: run_ctx

    ! scratch state vectors
    real(sp), allocatable :: state0(:), state1(:)

    ! differentiable work struct
    type(fuse_work) :: fuseStruct

  end type run_ctx

  CONTAINS
  
  SUBROUTINE fuse_evaluate(XPAR,GRID_FLAG,NCID_FORC,METRIC_VAL,OUTPUT_FLAG,IPSET,MPARAM_FLAG)

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
    use globaldata,only: NPAR_SNOW, isPrint, nFUSE_eval
    use model_defn,only: NSTATE
    use multiparam,only: NUMPAR
    use multiforce,only: nspat1, nspat2, numtim_sub
    use multibands,only: N_BANDS, n_bands
    use multistats,only: MSTATS, PCOUNT
    use multi_flux,only: W_FLUX_3d

    IMPLICIT NONE

    ! input
    REAL(SP),DIMENSION(:),INTENT(IN)       :: XPAR           ! model parameter set
    LOGICAL(LGT), INTENT(IN)               :: GRID_FLAG      ! .TRUE. if running FUSE on a grid
    INTEGER(I4B), INTENT(IN)               :: NCID_FORC      ! NetCDF ID for the forcing file
    LOGICAL(LGT), INTENT(IN)               :: OUTPUT_FLAG    ! .TRUE. if desire time series output
    INTEGER(I4B), INTENT(IN)               :: IPSET          ! index parameter set
    LOGICAL(LGT), INTENT(IN), OPTIONAL     :: MPARAM_FLAG    ! .FALSE. (used to turn off writing statistics)

    ! output
    REAL(SP),INTENT(OUT)                   :: METRIC_VAL     ! metric 

    ! run context
    type(run_ctx)                          :: ctx            ! container for allocatable structures

    ! error control
    integer(i4b)                       :: err, ierr
    character(len=1024)                :: message

    ! timing
    real(sp)                           :: t1, t2

    ! ---------------------------------------------------------------------------------------
    ! allocate run-time data structures
    call allocate_run(ctx, NSTATE, NUMPAR, N_BANDS, NPAR_SNOW, nspat1, nspat2, numtim_sub, ierr)
    if (ierr /= 0) stop "problem allocating run context in fuse_evaluate"

    ! allocate 3d data structure for fluxes
    allocate(w_flux_3d(nspat1, nspat2, numtim_sub), stat=ierr)
    if (ierr /= 0) stop "problem allocating w_flux_3d in fuse_evaluate"

    ! populate parameter structures and initialize states
    call initialize_run(ctx, XPAR, GRID_FLAG, MPARAM_FLAG, ierr, message)
    if (ierr /= 0) stop trim(message)
    
    ! initialize timing
    CALL CPU_TIME(T1)

    ! run fuse for the entire time series
    call run_time_loop(ctx, GRID_FLAG, NCID_FORC, OUTPUT_FLAG, err, message)
    if (err /= 0) stop trim(message)

    ! get timing information
    CALL CPU_TIME(T2)
    if(isPrint) WRITE(*,*) "TIME ELAPSED = ", t2-t1

    ! calculate mean summary statistics
    IF(.NOT.GRID_FLAG)THEN

      if(isPrint) PRINT *, 'Calculating performance metrics...'
      CALL MEAN_STATS()
      METRIC_VAL = MSTATS%METRIC_VAL

      write(*,'(i6,1x,a6,1x,f12.6,1x,a20,1x,f12.6)') nFUSE_eval, "NSE = ", MSTATS%NASH_SUTT, "; TIME ELAPSED = ", t2-t1
      !if(nFUSE_eval > 10) stop "checking results"

    ENDIF

    if(isPrint) PRINT *, 'Writing model statistics...'
    CALL PUT_SSTATS(PCOUNT)

    ! deallocate run context
    call deallocate_run(ctx, n_bands, ierr)
    if (ierr /= 0) stop "problem deallocating run context in fuse_evaluate"

    ! deallocate output buffer
    DEALLOCATE(W_FLUX_3d); IF (IERR.NE.0) STOP ' problem deallocating W_FLUX_3d in fuse_metric '

  END SUBROUTINE fuse_evaluate

  ! -------------------------------------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------------------------------------

  ! -------------------------------------------------------------------------------------------------------------------
  ! ----- private subroutine allocate_run: allocate run-time variables ------------------------------------------------
  ! -------------------------------------------------------------------------------------------------------------------

  subroutine allocate_run(ctx, nState, numpar, n_bands, npar_snow, nspat1, nspat2, numtim_sub, ierr)
  implicit none

  type(run_ctx), intent(inout) :: ctx
  integer(i4b), intent(in)     :: nState, numpar, n_bands, npar_snow, nspat1, nspat2, numtim_sub
  integer(i4b), intent(out)    :: ierr

  integer(i4b) :: iBands

  ierr = 0

  ! allocate state vectors
  allocate(ctx%state0(nState), ctx%state1(nState), stat=ierr)
  if (ierr /= 0) return

  ! allocate flux derivative vectors (inside fuseStruct)
  allocate(ctx%fuseStruct%df_dS(nState), ctx%fuseStruct%df_dPar(numpar), ctx%fuseStruct%dL_dPar(numpar), stat=ierr)
  if (ierr /= 0) return

  ! allocate elevation bands
  allocate(ctx%fuseStruct%sbands(n_bands), stat=ierr)
  if (ierr /= 0) return

  ! allocate parameter derivative for each elevation band
  do iBands = 1, n_bands

    allocate(ctx%fuseStruct%sbands(iBands)%var%dSWE_dParam(npar_snow), stat=ierr)
    if (ierr /= 0) return
  
    ctx%fuseStruct%sbands(iBands)%var%dSWE_dParam(:) = 0._sp
  
  end do

  end subroutine allocate_run

  ! -------------------------------------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------------------------------------

  ! -------------------------------------------------------------------------------------------------------------------
  ! ----- private subroutine deallocate_run: deallocate run-time variables --------------------------------------------
  ! -------------------------------------------------------------------------------------------------------------------

  subroutine deallocate_run(ctx, n_bands, ierr)
  implicit none

  type(run_ctx), intent(inout) :: ctx
  integer(i4b), intent(in)     :: n_bands
  integer(i4b), intent(out)    :: ierr
  integer(i4b)                 :: iBands

  ierr = 0

  ! deallocate parameter derivative vectors
  do iBands=1,n_bands
   deallocate(ctx%fuseStruct%sbands(iBands)%var%dSWE_dParam, stat=ierr)
   if (ierr /= 0) return
  end do

  ! deallocate state vectors
  DEALLOCATE(ctx%STATE0,ctx%STATE1,STAT=IERR)
  if (ierr /= 0) return
  
  ! deallocate flux derivative vectors
  deallocate(ctx%fuseStruct%df_dS, ctx%fuseStruct%df_dPar, ctx%fuseStruct%dL_dPar, stat=ierr)
  if (ierr /= 0) return

  ! deallocate elevation bands
  deallocate(ctx%fuseStruct%sbands, stat=ierr)
  if (ierr /= 0) return

  end subroutine deallocate_run

  ! -------------------------------------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------------------------------------

  ! -------------------------------------------------------------------------------------------------------------------
  ! ----- private subroutine initialize_run: populate param sets and initialize states  -------------------------------
  ! -------------------------------------------------------------------------------------------------------------------

  subroutine initialize_run(ctx, xpar, grid_flag, mparam_flag, err, message)
  
  use globaldata,  only: isPrint, fracstate0
  use model_defn,  only: SMODL
  use model_defnames
  
  use multiparam,  only: NUMPAR
  use multiforce,  only: nspat1, nspat2, DELTIM
  use multistate,  only: FSTATE, gState_3d
  use multistats,  only: PCOUNT
  use multibands
  
  use par_insert_module
  use str_2_xtry_module
  use xtry_2_str_module
  use put_params_module, only: put_params
  implicit none

  type(run_ctx), intent(inout)               :: ctx
  real(sp), dimension(:), intent(in)         :: xpar
  logical(lgt), intent(in)                   :: grid_flag
  logical(lgt), intent(in), optional         :: mparam_flag

  integer(i4b), intent(out)                  :: err
  character(len=*), intent(out)              :: message

  integer(i4b) :: iSpat1, iSpat2, iBands

  err = 0
  message = ""

  ! increment parameter counter for model output
  if (.not. present(mparam_flag)) then
    PCOUNT = PCOUNT + 1
  else
    if (mparam_flag) PCOUNT = PCOUNT + 1
  end if

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
  if (SMODL%iSNOWM == iopt_temp_index .and. .not. grid_flag) then
    Z_FORCING      = Z_FORCING_grid(1,1)
    MBANDS(:)%info = MBANDS_INFO_3d(1,1,:)
  end if

  if (isPrint) print *, 'Writing parameter values...'
  call put_params(PCOUNT)

  ! initialize model states over the 2D gridded domain (1x1 in catchment mode)
  do iSpat2 = 1, nSpat2
    do iSpat1 = 1, nSpat1
      call init_state(fracstate0)
      call str_2_xtry(FSTATE, ctx%state0)
      call xtry_2_str(ctx%state0, FSTATE)
      gState_3d(iSpat1, iSpat2, 1) = FSTATE
    end do
  end do
  if (isPrint) print *, 'Model states initialized over the 2D gridded domain'

  ! initialize elevation bands if snow module is on
  if (isPrint) print *, 'N_BANDS =', N_BANDS
  if (SMODL%iSNOWM == iopt_temp_index) then

    ! initialize template once
    ctx%fuseStruct%sbands(:)%var%SWE         = 0._sp
    ctx%fuseStruct%sbands(:)%var%SNOWACCMLTN = 0._sp
    ctx%fuseStruct%sbands(:)%var%SNOWMELT    = 0._sp
    ctx%fuseStruct%sbands(:)%var%DSWE_DT     = 0._sp

    ! copy to every grid cell (legacy staging)
    do iSpat2 = 1, nSpat2
      do iSpat1 = 1, nSpat1
        do iBands = 1, n_bands
          MBANDS_VAR_4d(iSpat1, iSpat2, iBands, 1) = ctx%fuseStruct%sbands(iBands)%var%bands_var
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

  subroutine run_time_loop(ctx, grid_flag, ncid_forc, output_flag, err, message)

  ! switches / options
  use globaldata, only: isPrint, NA_VALUE_SP
  use model_defn,  only: SMODL, NSTATE
  use model_defnames
  use multiforce, only: &
      nspat1, nspat2, DELTIM, &
      sim_beg, sim_end, &
      numtim_sub, numtim_sim, numtim_sub_cur, &
      itim_in, itim_sim, itim_sub, &
      gForce_3d, aForce, MFORCE
  use multistate, only: gState_3d, FSTATE, MSTATE
  use multiroute, only: MROUTE, AROUTE_3d
  use multibands
  use multi_flux, only: W_FLUX, W_FLUX_3d
  
  ! modules
  use time_io,           only: get_modtim
  use get_gforce_module, only: get_gforce_3d
  use getPETgrid_module, only: getPETgrid
  use put_output_module, only: put_goutput_3d
  use set_all_module,    only: SET_STATE, SET_FLUXES, SET_ROUTE
  
  ! state vector conversions
  use str_2_xtry_module
  use xtry_2_str_module
  
  ! differentiable
  use get_bundle_module,      only: get_bundle
  use implicit_solve_module,  only: implicit_solve
  use update_swe_diff_module, only: update_swe_diff   ! confirm spelling matches your subroutine

  ! original solver interface
  use interfaceb, only: ode_int, fuse_solve

  ! model numerix structures
  USE model_numerix
  USE fuse_deriv_module
  USE fdjac_ode_module

  implicit none

  type(run_ctx), intent(inout) :: ctx
  logical(lgt),  intent(in)    :: grid_flag
  integer(i4b),  intent(in)    :: ncid_forc
  logical(lgt),  intent(in)    :: output_flag

  integer(i4b),  intent(out)   :: err
  character(len=*), intent(out):: message

  ! locals
  logical(lgt), parameter :: computePET = .false.
  real(sp)     :: dt_sub, dt_full
  integer(i4b) :: iSpat1, iSpat2, iBands
  integer(i4b) :: ierr
  character(len=1024) :: cmessage

  err = 0
  message = ""
  cmessage = ""
  ierr = 0

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
  itim_sub = 1
  itim_sim = 1

  ! loop through time steps of the input file (ITIM_IN)
  DO ITIM_IN=sim_beg,sim_end

    ! if start of subperiod: load forcing
    IF(itim_sub.EQ.1)THEN

      ! determine length of current subperiod
      numtim_sub_cur=MIN(numtim_sub,numtim_sim-itim_sim+1)

      ! load forcing for desired period into gForce_3d
      if(isPrint) PRINT *, 'New subperiod: loading forcing for ',numtim_sub_cur,' time steps'
      CALL get_gforce_3d(itim_in,numtim_sub_cur,ncid_forc,err,message)
      IF(err/=0)THEN; WRITE(*,*) 'Error while extracting 3d forcing'; STOP; ENDIF
      if(isPrint) PRINT *, 'Forcing loaded. Running FUSE...'

    ENDIF

    ! get the model time
    CALL get_modtim(itim_in,ncid_forc,ierr,message)
    IF(ierr/=0)THEN; PRINT*, TRIM(cmessage); STOP; ENDIF
    !print*, timdat

    ! compute potential ET
    IF(computePET) CALL getPETgrid(ierr,cmessage)
    IF(ierr/=0)THEN; PRINT*, TRIM(cmessage); STOP; ENDIF

    ! loop through grid points and run the model for one time step
    DO iSpat2=1,nSpat2
      DO iSpat1=1,nSpat1

          ! only run FUSE for grid points within domain defined by elev_mask
          IF(.NOT.elev_mask(iSpat1,iSpat2))THEN

            ! FUSE works with MFORCE, MSTATE, MBANDS, W_FLUX, MROUTE, which are all scalars.
            ! Here we transfer forcing, state, flux variables from the 3D structures to these
            ! variables, run FUSE and then transfer the new values back to the 3D structures.

            ! extract forcing for this grid cell and time step
            MFORCE = gForce_3d(iSpat1,iSpat2,itim_sub)

            ! forcing sanity checks
            if(MFORCE%PPT.lt.0.0) then; PRINT *, 'Negative precipitation in input file:',iSpat1,iSpat2,MFORCE%PPT; stop; endif
            if(MFORCE%PPT.gt.5000.0) then; PRINT *, 'Precipitation greater than 5000 in input file:',iSpat1,iSpat2,MFORCE%PPT; stop; endif
            if(MFORCE%PET.lt.0.0) then; PRINT *, 'Negative PET in input file'; stop; endif
            if(MFORCE%PET.gt.100.0) then; PRINT *, 'PET greater than 100 in input file'; stop; endif
            if(MFORCE%TEMP.lt.-100.0) then; PRINT *, 'Temperature lower than -100 in input file'; stop; endif
            if(MFORCE%TEMP.gt.100.0) then; PRINT *, 'Temperature greater than 100 in input file'; stop; endif

             ! extract model states for this grid cell and time step
             FSTATE = gState_3d(iSpat1,iSpat2,itim_sub)
             MSTATE = FSTATE                     ! refresh model states
             CALL STR_2_XTRY(FSTATE, ctx%STATE0) ! set state at the start of the time step (STATE0) using FSTATE

             ! initialize model fluxes
             CALL INITFLUXES()                   ! set weighted sum of fluxes to zero

             ! populate fuse work structure
             if(diff_mode==differentiable) call get_bundle(ctx%fuseStruct)

             ! if snow model is on, call UPDATE_SWE to calculate snow fluxes and update snow bands
             ! using explicit Euler approach; if not, call QRAINERROR
             SELECT CASE(SMODL%iSNOWM)
             CASE(iopt_temp_index)

                ! load data from multidimensional arrays
                Z_FORCING      = Z_FORCING_grid(iSpat1,iSpat2)             ! elevation of forcing data (m)
                mbands(:)%info = MBANDS_INFO_3d(iSpat1,iSpat2,:)           ! info structure
                mbands(:)%var  = MBANDS_VAR_4d(iSpat1,iSpat2,:,itim_sub)   ! var structure

                ! put data into the FUSE structure
                ! NOTE: only copy the "var" variables
                if(diff_mode == differentiable)then
                 ctx%fuseStruct%z_forcing                = Z_FORCING
                 ctx%fuseStruct%sbands(:)%info           = MBANDS(:)%info
                 ctx%fuseStruct%sbands(:)%var%bands_var  = MBANDS(:)%var
                endif  ! if diff_mode == differentiable

                ! run the snow model
                select case(diff_mode)
                 case(original);       CALL UPDATE_SWE(DELTIM)
                 case(differentiable); CALL UPDATE_SWE_DIFF(ctx%fuseStruct,DELTIM)
                 CASE DEFAULT; stop "fuse_metric: cannot identify diff_mode"
                end select

             CASE(iopt_no_snowmod)
                CALL QRAINERROR()
             CASE DEFAULT
                message="fuse_metric/SMODL%iSNOWM must be either iopt_temp_index or iopt_no_snowmod"
                print*, trim(message); stop 1
             END SELECT

            ! ----- start of soil physics code ------------------------------------------------------------

            ! temporally integrate the ordinary differential equations
            select case(diff_mode)

              ! original code
              case(original)
               CALL ODE_INT(FUSE_SOLVE,ctx%STATE0,ctx%STATE1,DT_SUB,DT_FULL,IERR,MESSAGE)
               IF (IERR.NE.0) THEN; PRINT *, TRIM(MESSAGE); STOP 1; ENDIF

              ! differentiable code
              case(differentiable)

               ! solve differentiable ODEs
               call implicit_solve(ctx%fuseStruct, ctx%state0, ctx%state1, nState, ierr, cmessage)

               ! save fluxes
               W_FLUX = ctx%fuseStruct%flux

              ! check options
              case default; print*, "fuse_metric: Cannot identify diff_mode"; stop 1
            end select

            !print*, ITIM_IN, w_flux%eff_ppt
            !if(ITIM_IN > 100) stop "check"

            ! ----- end of soil physics code --------------------------------------------------------------

            ! perform overland flow routing
            CALL Q_OVERLAND()

            ! runoff sanity check
            IF (MROUTE%Q_ROUTED.LT.0._sp) STOP 'Q_ROUTED is less than zero'
            IF (MROUTE%Q_ROUTED.GT.1000._sp) STOP 'Q_ROUTED is enormous'

            ! transfer simulations to corresponding 3D structures
            ! note that the first time step of gState_3d and MBANDS_VAR_4d is defined by initialisation
            ! or simulation over previous subperiod, so saving in itim_sub+1 - and hence, the allocated
            ! length of the temporal dimension of gState_3d and MBANDS_VAR_4d is numtim_sub+1,
            ! but numtim_sub for W_FLUX_3d and AROUTE_3d

            CALL XTRY_2_STR(ctx%STATE1,FSTATE)            ! update FSTATE using states at the end of the time step (STATE1)
            gState_3d(iSpat1,iSpat2,itim_sub+1) = FSTATE  ! transfer FSTATE into the 3-d structure
            W_FLUX_3d(iSpat1,iSpat2,itim_sub) = W_FLUX    ! fluxes
            AROUTE_3d(iSpat1,iSpat2,itim_sub) = MROUTE    ! instantaneous and routed runoff

            IF (SMODL%iSNOWM.EQ.iopt_temp_index) THEN

              ! extract data from the FUSE structure
              if(diff_mode == differentiable)then
               Z_FORCING   = ctx%fuseStruct%z_forcing
               MBANDS%info = ctx%fuseStruct%sbands%info         
               MBANDS%var  = ctx%fuseStruct%sbands%var%bands_var
              endif  ! if diff_mode == differentiable

              ! SWE TOT: weighted average of SWE over all the elevation bands
              gState_3d(iSpat1,iSpat2,itim_sub+1)%SWE_TOT = SUM(MBANDS(:)%var%SWE * MBANDS(:)%info%AF)

              ! update MBANDS_VAR_4D
              MBANDS_VAR_4d(iSpat1,iSpat2,:,itim_sub+1) = MBANDS(:)%var

            END IF

            ! save forcing data to export to output file
            IF(GRID_FLAG)THEN
               aForce(itim_sub)%ppt = SUM(gForce_3d(:,:,itim_sub)%ppt)/REAL(SIZE(gForce_3d(:,:,itim_sub)), KIND(sp))
               aForce(itim_sub)%pet = SUM(gForce_3d(:,:,itim_sub)%pet)/REAL(SIZE(gForce_3d(:,:,itim_sub)), KIND(sp))
            ENDIF

            ! compute summary statistics
            CALL COMP_STATS()

          ELSE ! insert NA values if grid point outside of domain or forcing not available

            CALL SET_STATE(NA_VALUE_SP) ! includes FSTATE%SWE_TOT
            gState_3d(iSpat1,iSpat2,itim_sub) = FSTATE

            CALL SET_FLUXES(NA_VALUE_SP)
            W_FLUX_3d(iSpat1,iSpat2,itim_sub) = W_FLUX

            CALL SET_ROUTE(NA_VALUE_SP)
            AROUTE_3d(iSpat1,iSpat2,itim_sub) = MROUTE

         ENDIF ! (is grid cell in mask_elev?)
      END DO  ! (looping thru 2nd spatial dimension)
    END DO  ! (looping thru 1st spatial dimension)

    ! if end of subperiod: write to output file and save states
    IF(itim_sub.EQ.numtim_sub_cur)THEN

      if(isPrint) PRINT *, 'End of subperiod reached:'

      ! write model output
      IF (OUTPUT_FLAG) THEN
        if(isPrint) PRINT *, 'Write output for ',numtim_sub_cur,' time steps starting at indices', itim_sim-numtim_sub_cur+1
        CALL PUT_GOUTPUT_3D(itim_sim-numtim_sub_cur+1, itim_in-numtim_sub_cur+1, numtim_sub_cur)
        if(isPrint) PRINT *, 'Done writing output'
      ELSE
        if(isPrint) PRINT *, 'OUTPUT_FLAG is set on FALSE, no output written'
      END IF

      ! TODO: set gState_3d and MBANDS_VAR_4d to NA

      ! reinitialize states for next subperiod using last time step
      gState_3d(:,:,1)       = gState_3d(:,:,itim_sub+1)
      MBANDS_VAR_4d(:,:,:,1) = MBANDS_VAR_4d(:,:,:,itim_sub+1)

      ! reset itim_sub
      itim_sub=1

    ELSE ! not the end of subperiod

      ! increment itim_sub
      itim_sub=itim_sub+1

    END IF

  ! increment itim_sim
  itim_sim=itim_sim+1

  END DO  ! (loop through timesteps)

  end subroutine run_time_loop
  
  ! -------------------------------------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------------------------------------

  ! -------------------------------------------------------------------------------------------------------------------
  ! ----- private subroutine 
  ! -------------------------------------------------------------------------------------------------------------------












END MODULE fuse_evaluate_module
