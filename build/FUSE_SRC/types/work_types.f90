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
 public :: fuse_work

 ! --------------------------------------------------------------------------------------
 
 ! dSWE/dParam for each elevation band

 type, extends(bands_var) :: bands_var_diff
   real(sp), allocatable               :: dSWE_dParam(:)
 end type bands_var_diff

 ! extended bands structure
 type ebands
   type(bands_info)                    :: info
   type(bands_var_diff)                :: var
 end type ebands

 ! --------------------------------------------------------------------------------------
 
 ! omnibus structure that bundles "everything" required to run fuse for a single cell

 type fuse_work
   type(tdata)                         :: time           ! time data
   type(fdata)                         :: force          ! model forcing data
   type(ebands) , allocatable          :: sbands(:)      ! info/variables for elevation bands (snow model)
   type(statev)                        :: state0         ! state variables (start of step)
   type(statev)                        :: state1         ! state variables (end of step)
   type(statev)                        :: dx_dt          ! time derivative in state variables
   type(fluxes)                        :: flux           ! fluxes
   type(fluxes), allocatable           :: df_dS(:)       ! derivative in fluxes w.r.t. states
   type(fluxes), allocatable           :: df_dPar(:)     ! derivative in fluxes w.r.t. parameters
   real(sp),     allocatable           :: dL_dPar(:)     ! derivative in loss function w.r.t. parameters
   type(runoff)                        :: route          ! hillslope routing
   type(par_id)                        :: param_name     ! parameter names
   type(parinfo)                       :: param_meta     ! metadata on model parameters
   type(paradj)                        :: param_adjust   ! adjustable model parametrs
   type(pardvd)                        :: param_derive   ! derived model parameters
   type(summary)                       :: sim_stats      ! simulation statistics
   real(sp)                            :: z_forcing      ! elevation of forcing data (m)
   logical(lgt)                        :: is_initialized = .false.
 end type fuse_work

end module work_types
