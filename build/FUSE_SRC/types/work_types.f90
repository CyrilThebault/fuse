module work_types

 ! data types

 use nrtype

 use multiforce_types,  only: TDATA, VDATA, ADATA, FDATA
 use multibands_types,  only: BANDS, BANDS_INFO, BANDS_VAR
 use multiparam_types,  only: PARATT, PARINFO, PARADJ, PARDVD, PAR_ID
 use multistate_types,  only: STATEV, M_TIME
 use multi_flux_types,  only: FLUXES
 use multiroute_types,  only: RUNOFF

 use multistats_types,  only: SUMMARY

 private

 public :: bands_var_diff, ebands
 public :: fuse_param
 public :: fuse_work

 ! --------------------------------------------------------------------------------------
 
 ! dSWE/dParam for each elevation band

 type, extends(bands_var) :: bands_var_diff
   real(wp), allocatable               :: dSWE_dParam(:)
 end type bands_var_diff

 ! extended bands structure
 type ebands
   type(bands_info)                    :: info
   type(bands_var_diff)                :: var
 end type ebands

 ! --------------------------------------------------------------------------------------
 ! structure bundles

 ! per-step structure
 type fuse_step
   type(tdata)                         :: time           ! time data
   type(fdata)                         :: force          ! model forcing data
   type(statev)                        :: state0         ! state variables (start of step)
   type(statev)                        :: state1         ! state variables (end of step)
   type(statev)                        :: tState         ! state variables (trial)
   type(statev)                        :: dx_dt          ! time derivative in state variables
   type(fluxes)                        :: flux           ! fluxes
   type(runoff)                        :: route          ! hillslope routing
 end type fuse_step

 ! snow structure
 type fuse_snow
   real(wp)                            :: z_forcing      ! elevation of forcing data (m)
   type(ebands) , allocatable          :: sbands(:)      ! info/variables for elevation bands (snow model)
 end type fuse_snow

 ! parameter structure
 type fuse_param
   type(par_id)                        :: param_name     ! parameter names
   type(parinfo)                       :: param_meta     ! metadata on model parameters
   type(paradj)                        :: param_adjust   ! adjustable model parametrs
   type(pardvd)                        :: param_derive   ! derived model parameters
 end type fuse_param

 ! numerix structure (linear algebra, ...)
 type fuse_numerix
   real(wp)     , allocatable          :: x0(:)      ! state variables (start of step)
   real(wp)     , allocatable          :: x1(:)      ! state variables (end of step)
 end type fuse_numerix

 ! adjoint structure (differentiable fuse)
 type fuse_adjoint
   type(fluxes), allocatable           :: df_dS(:)       ! derivative in fluxes w.r.t. states
   type(fluxes), allocatable           :: df_dPar(:)     ! derivative in fluxes w.r.t. parameters
   real(wp),     allocatable           :: dL_dPar(:)     ! derivative in loss function w.r.t. parameters
 end type fuse_adjoint

 ! run-level / evaluation-level
 type fuse_run
   type(summary) :: stats
   integer(i4b)  :: n_evaluations = 0  ! number of parameter sets evaluated
 end type fuse_run

 ! --------------------------------------------------------------------------------------
 ! omnibus structure that bundles "everything" required to run fuse for a single cell
 
 type fuse_work
   type(fuse_step)    :: step    ! per-step structure
   type(fuse_snow)    :: snow    ! snow structure
   type(fuse_param)   :: par     ! parameter structure
   type(fuse_numerix) :: num     ! numerix structure (linear algebra, ...)
   type(fuse_adjoint) :: adj     ! adjoint structure (differentiable fuse)
   type(fuse_run)     :: run     ! run-level structure
   logical(lgt)       :: is_initialized = .false.
 end type fuse_work

end module work_types
