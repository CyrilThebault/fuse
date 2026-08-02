! ---------------------------------------------------------------------
! Compatibility modules for compiling selected mizuRoute routines in FUSE.
!
! These modules intentionally preserve the module names expected by
! unmodified mizuRoute source files. They expose only the minimal subset
! of symbols required by the FUSE integration and are not complete
! replacements for the corresponding mizuRoute modules.
!
! Do not add new FUSE functionality here. This file is transitional and
! should remain as small as possible.
! ---------------------------------------------------------------------

!-----------------------------------------------------------------------
! Compatibility shim for mizuRoute.
! Provides the NetCDF data-type identifiers required by popMetadat
! without introducing a dependency on the PIO library.
!-----------------------------------------------------------------------
module pio_utils

  use nrtype, only: i4b
  use netcdf, only: nf90_float, nf90_double, nf90_int

  implicit none
  private

  public :: ncd_float
  public :: ncd_double
  public :: ncd_int

  integer(i4b), parameter :: ncd_float  = nf90_float
  integer(i4b), parameter :: ncd_double = nf90_double
  integer(i4b), parameter :: ncd_int    = nf90_int

end module pio_utils

!-----------------------------------------------------------------------
! Compatibility shim for mizuRoute.
! Provides the minimal subset of global variables required by the
! unmodified mizuRoute topology routines compiled within FUSE.
! These include metadata describing the river-network data structures
! and routing-method flags used during topology initialization.
!-----------------------------------------------------------------------
module globalData

  ! Thin compatibility module for the initial mizuRoute foundation.
  ! Contains only the global metadata and process flag required by
  ! read_streamSeg and mod_meta_varFile.

  use nrtype,     only: i4b, dp, lgt
  
  USE dataTypes, ONLY: struct_info   ! metadata type - data structure
  USE dataTypes, ONLY: dim_info      ! metadata type - variable dimensions
  USE dataTypes, ONLY: var_info      ! metadata type - variable
  USE objTypes,  ONLY: var_info_new  ! metadata type - variable

  USE var_lookup, ONLY: nStructures  ! number of variables in data structure (struct_info)
  USE var_lookup, ONLY: nDimensions  ! number of variables in dimensions related to network topology
  USE var_lookup, ONLY: nStateDims   ! number of variables in dimensions related to restart variables
  USE var_lookup, ONLY: nQdims       ! number of variables in dimensions related to fluxes/states variables
  USE var_lookup, ONLY: nVarsHRU     ! number of variables in data structure (catchment propoerties)
  USE var_lookup, ONLY: nVarsHRU2SEG ! number of variables in data structure (river-catchment topology)
  USE var_lookup, ONLY: nVarsSEG     ! number of variables in data structure (river reach propeties)
  USE var_lookup, ONLY: nVarsNTOPO   ! number of variables in data structure (river network topology)
  USE var_lookup, ONLY: nVarsPFAF    ! number of variables in data structure (pfaffstetter related variable)
  USE var_lookup, ONLY: nVarsRFLX    ! number of variables in data structure (river flux/state)
  USE var_lookup, ONLY: nVarsHFLX    ! number of variables in data structure (HRU flux/state)
  USE var_lookup, ONLY: nVarsBasinQ  ! number of variables in data structure (restart vars for
  USE var_lookup, ONLY: nVarsIRFbas  ! number of variables in data structure (restart vars for overland unit-hydrograph routing)
  USE var_lookup, ONLY: nVarsBasTracer ! number of variables in data structure (restart vars for overland tracer routing)
  USE var_lookup, ONLY: nVarsIRF     ! number of variables in data structure (restart vars for unit-hydrograph routing)
  USE var_lookup, ONLY: nVarsKWT     ! number of variables in data structure (restart vars for lagrangian kinematic wave)
  USE var_lookup, ONLY: nVarsKW      ! number of variables in data structure (restart vars for kinematic wave routing)
  USE var_lookup, ONLY: nVarsMC      ! number of variables in data structure (restart vars for muskingum-cunge routing)
  USE var_lookup, ONLY: nVarsDW      ! number of variables in data structure (restart vars for diffusive wave routing)
  USE var_lookup, ONLY: nVarsTracer  ! number of variables in data structure (restart vars for tracer)

  use public_var, only: nRouteMethods
  
  implicit none
  private
  save

  ! ---------- Misc. data -------------------------------------------------------------------------

  logical(lgt),                    public :: masterproc
  integer(i4b),                    public :: maxtdh=0                    ! maximum unit-hydrograph future time steps

  ! time delay histogram (hillslope routing)
  real(dp),           allocatable, public :: FRAC_FUTURE(:)         ! fraction of runoff in future time steps

  ! logical to indicate active routing method(s)
  !  -- named indices for routing methods defined in public_var
  logical(lgt), public :: onRoute(0:nRouteMethods-1) = .false.

  ! ---------- routing parameter names -------------------------------------------------------------------
  ! spatially constant ....
  real(dp),                        public :: fshape                     ! shape parameter in time delay histogram (=gamma distribution) [-]
  real(dp),                        public :: tscale                     ! scaling factor for the time delay histogram [sec]
  real(dp),                        public :: velo                       ! velocity [m/s] for Saint-Venant equation
  real(dp),                        public :: diff                       ! diffusivity [m2/s] for Saint-Venant equation
  real(dp),                        public :: mann_n                     ! manning's roughness coefficient [-]
  real(dp),                        public :: wscale                     ! scaling factor for river width [-]
  real(dp),                        public :: dscale=0.000045            ! scaling factor for river bankful depth [-]
  real(dp),                        public :: floodplainSlope=1000       ! floodplain down slope h:v=slope:1 [-]
  real(dp),                        public :: high_depth=100000._dp      ! very high river bankful depth [m]

  ! ---------- general structure information --------------------------------------------------------

  type(struct_info),               public :: meta_struct(nStructures)   ! metadata on the data structures
  type(dim_info),                  public :: meta_dims(nDimensions)     ! metadata on the dimensions for network topology
  type(dim_info),                  public :: meta_stateDims(nStateDims) ! metadata on the dimensions for state variables
  type(dim_info),                  public :: meta_qDims(nQdims)         ! metadata on the dimensions for flux variables
  type(dim_info),                  public :: meta_qDims_gage(nQdims)    ! metadata on the dimensions for flux variables


  ! ---------- metadata structures ------------------------------------------------------------------
  
  type(var_info),                  public :: meta_HRU    (nVarsHRU    )      ! HRU properties
  type(var_info),                  public :: meta_HRU2SEG(nVarsHRU2SEG)      ! HRU-to-segment mapping
  type(var_info),                  public :: meta_SEG    (nVarsSEG    )      ! stream segment properties
  type(var_info),                  public :: meta_NTOPO  (nVarsNTOPO  )      ! network topology
  type(var_info),                  public :: meta_PFAF   (nVarsPFAF   )      ! pfafstetter code
  type(var_info_new),              public :: meta_rflx   (nVarsRFLX   )      ! reach flux variables
  type(var_info_new),              public :: meta_hflx   (nVarsHFLX   )      ! hru flux variables
  type(var_info_new),              public :: meta_basinQ (nVarsBasinQ )      ! reach inflow from basin
  type(var_info_new),              public :: meta_irf_bas(nVarsIRFbas )      ! basin IRF routing fluxes/states
  type(var_info_new),              public :: meta_bas_solute(nVarsBasTracer) ! basin IRF routing for solute mass flux
  type(var_info_new),              public :: meta_irf    (nVarsIRF    )      ! IRF routing fluxes/states
  type(var_info_new),              public :: meta_kwt    (nVarsKWT    )      ! KWT routing fluxes/states
  type(var_info_new),              public :: meta_kw     (nVarsKW     )      ! KW routing fluxes/states
  type(var_info_new),              public :: meta_mc     (nVarsMC     )      ! MC routing restart fluxes/states
  type(var_info_new),              public :: meta_dw     (nVarsDW     )      ! DW routing restart fluxes/states
  type(var_info_new),              public :: meta_solute (nVarsTracer )      ! solute mass fluxes states

end module globalData
