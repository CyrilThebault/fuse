MODULE multiparam_types
 
 ! ---------------------------------------------------------------------------------------
 ! Creator:
 ! --------
 ! Martyn Clark
 ! Modified by Brian Henn to include snow model, 6/2013
 ! Modified by Martyn Clark to separate type definitions from data storage, 01/2026
 ! ---------------------------------------------------------------------------------------
 
 USE nrtype
 USE model_defn, ONLY: NTDH_MAX
 
 implicit none
 private

 public :: PARATT, PARINFO, PARADJ, PARDVD, PAR_ID

 ! --------------------------------------------------------------------------------------
 ! (1) PARAMETER METADATA
 ! --------------------------------------------------------------------------------------
 ! data structure to hold metadata for adjustable model parameters
 TYPE PARATT
  LOGICAL(LGT)                         :: PARFIT      ! flag to determine if parameter is fitted
  INTEGER(I4B)                         :: PARSTK      ! flag (0=deterministic, 1=stochastic)
  REAL(WP)                             :: PARDEF      ! default parameter set
  REAL(WP)                             :: PARLOW      ! lower limit of each parameter
  REAL(WP)                             :: PARUPP      ! upper limit of each parameter
  REAL(WP)                             :: FRSEED      ! fraction param space for "reasonable" bounds
  REAL(WP)                             :: PARSCL      ! typical scale of parameter
  INTEGER(I4B)                         :: PARVTN      ! method used for variable transformation
  INTEGER(I4B)                         :: PARDIS      ! parametric form of prob dist used for prior/hyper
  INTEGER(I4B)                         :: PARQTN      ! transformation applied before use of prob dist
  INTEGER(I4B)                         :: PARLAT      ! number of latent variables (0=onePerStep, -1=from data)
  INTEGER(I4B)                         :: PARMTH      ! imeth for all variables ???what is this???
  INTEGER(I4B)                         :: NPRIOR      ! number of prior/hyper-parameters
  CHARACTER(LEN=256)                   :: P_NAME      ! parameter name
  CHARACTER(LEN=256)                   :: CHILD1      ! name of 1st parameter child
  CHARACTER(LEN=256)                   :: CHILD2      ! name of 2nd parameter child
 END TYPE PARATT
 
 ! data structure to hold metadata for each parameter
 TYPE PARINFO
  ! rainfall error parameters (adjustable)
  TYPE(PARATT)                         :: RFERR_ADD   !  additive rainfall error (mm day-1)
  TYPE(PARATT)                         :: RFERR_MLT   ! multiplicative rainfall error (-)
  TYPE(PARATT)                         :: RFH1_MEAN   ! hyper parameter1: mean rainfall multiplier (-)
  TYPE(PARATT)                         :: RFH2_SDEV   ! hyper parameter2: sdev rainfall multiplier (-)
  TYPE(PARATT)                         :: RH1P_MEAN   ! prior param1 of hyper param1: prior mean of hypermean
  TYPE(PARATT)                         :: RH1P_SDEV   ! prior param2 of hyper param1: prior sdev of hypermean
  TYPE(PARATT)                         :: RH2P_MEAN   ! prior param1 of hyper param2: lower bound of hypersdev
  TYPE(PARATT)                         :: RH2P_SDEV   ! prior param2 of hyper param2: upper bound of hypersdev
  ! bucket sizes (adjustable)
  TYPE(PARATT)                         :: MAXWATR_1   ! maximum total storage in layer1 (mm)
  TYPE(PARATT)                         :: MAXWATR_2   ! maximum total storage in layer2 (mm)
  TYPE(PARATT)                         :: FRACTEN     ! frac total storage as tension storage (-)
  TYPE(PARATT)                         :: FRCHZNE     ! PRMS: frac tension storage in recharge zone (-)
  TYPE(PARATT)                         :: FPRIMQB     ! SAC: fraction of baseflow in primary resvr (-)
  ! evaporation (adjustable)
  TYPE(PARATT)                         :: RTFRAC1     ! fraction of roots in the upper layer (-)
  ! percolation (adjustable) 
  TYPE(PARATT)                         :: PERCRTE     ! percolation rate (mm day-1)
  TYPE(PARATT)                         :: PERCEXP     ! percolation exponent (-)
  TYPE(PARATT)                         :: SACPMLT     ! multiplier in the SAC model for dry lower layer (-)
  TYPE(PARATT)                         :: SACPEXP     ! exponent in the SAC model for dry lower layer (-)
  TYPE(PARATT)                         :: PERCFRAC    ! fraction of percolation to tension storage (-)
  TYPE(PARATT)                         :: FRACLOWZ    ! fraction of soil excess to lower zone (-)
  ! interflow (adjustable)
  TYPE(PARATT)                         :: IFLWRTE     ! interflow rate (mm day-1)
  ! baseflow (adjustable)
  TYPE(PARATT)                         :: BASERTE     ! baseflow rate (mm day-1)
  TYPE(PARATT)                         :: QB_POWR     ! baseflow exponent (-)
  TYPE(PARATT)                         :: QB_PRMS     ! baseflow depletion rate (day-1)
  TYPE(PARATT)                         :: QBRATE_2A   ! baseflow depletion rate for primary resvr (day-1)
  TYPE(PARATT)                         :: QBRATE_2B   ! baseflow depletion rate for secondary resvr (day-1)
  ! surface runoff (adjustable)
  TYPE(PARATT)                         :: SAREAMAX    ! maximum saturated area
  TYPE(PARATT)                         :: AXV_BEXP    ! ARNO/VIC "b" exponent
  TYPE(PARATT)                         :: LOGLAMB     ! mean value of the log-transformed topographic index (m)
  TYPE(PARATT)                         :: TISHAPE     ! shape parameter for the topo index Gamma distribution (-)
  ! time delay in runoff
  TYPE(PARATT)                         :: TIMEDELAY   ! time delay in runoff (days)
  ! snow model (adjustable)
  TYPE(PARATT)                         :: MBASE       ! base melt temperature (deg. C)
  TYPE(PARATT)                         :: MFMAX       ! maximum melt factor (mm melt deg C.-1 6hrs-1)
  TYPE(PARATT)                         :: MFMIN       ! minimum melt factor (mm melt deg C.-1 6hrs-1)
  TYPE(PARATT)                         :: PXTEMP      ! rain-snow partition temperature (deg. C)
  TYPE(PARATT)                         :: OPG         ! precipitation gradient (-) 
  TYPE(PARATT)                         :: LAPSE       ! temperature gradient (deg. C)
 ENDTYPE PARINFO
 
 ! --------------------------------------------------------------------------------------
 ! (2) ADJUSTABLE PARAMETERS
 ! --------------------------------------------------------------------------------------
 TYPE PARADJ
  ! rainfall error parameters (adjustable)
  REAL(WP)                             :: RFERR_ADD   ! additive rainfall error (mm day-1)
  REAL(WP)                             :: RFERR_MLT   ! multiplicative rainfall error (-)
  REAL(WP)                             :: RFH1_MEAN   ! hyper parameter1: mean rainfall multiplier (-)
  REAL(WP)                             :: RFH2_SDEV   ! hyper parameter2: sdev rainfall multiplier (-)
  REAL(WP)                             :: RH1P_MEAN   ! prior param1 of hyper param1: prior mean of hypermean
  REAL(WP)                             :: RH1P_SDEV   ! prior param2 of hyper param1: prior sdev of hypermean
  REAL(WP)                             :: RH2P_MEAN   ! prior param1 of hyper param2: lower bound of hypersdev
  REAL(WP)                             :: RH2P_SDEV   ! prior param2 of hyper param2: upper bound of hypersdev
  ! bucket sizes (adjustable)
  REAL(WP)                             :: MAXWATR_1   ! maximum total storage in layer1 (mm)
  REAL(WP)                             :: MAXWATR_2   ! maximum total storage in layer2 (mm)
  REAL(WP)                             :: FRACTEN     ! frac total storage as tension storage (-)
  REAL(WP)                             :: FRCHZNE     ! PRMS: frac tension storage in recharge zone (-)
  REAL(WP)                             :: FPRIMQB     ! SAC: fraction of baseflow in primary resvr (-)
  ! evaporation (adjustable)
  REAL(WP)                             :: RTFRAC1     ! fraction of roots in the upper layer (-)
  ! percolation (adjustable)
  REAL(WP)                             :: PERCRTE     ! percolation rate (mm day-1)
  REAL(WP)                             :: PERCEXP     ! percolation exponent (-)
  REAL(WP)                             :: SACPMLT     ! multiplier in the SAC model for dry lower layer (-)
  REAL(WP)                             :: SACPEXP     ! exponent in the SAC model for dry lower layer (-)
  REAL(WP)                             :: PERCFRAC    ! fraction of percolation to tension storage (-)
  REAL(WP)                             :: FRACLOWZ    ! fraction of soil excess to lower zone (-)
  ! interflow (adjustable)
  REAL(WP)                             :: IFLWRTE     ! interflow rate (mm day-1)
  ! baseflow (adjustable)
  REAL(WP)                             :: BASERTE     ! baseflow rate (mm day-1)
  REAL(WP)                             :: QB_POWR     ! baseflow exponent (-)
  REAL(WP)                             :: QB_PRMS     ! baseflow depletion rate (day-1)
  REAL(WP)                             :: QBRATE_2A   ! baseflow depletion rate for primary resvr (day-1)
  REAL(WP)                             :: QBRATE_2B   ! baseflow depletion rate for secondary resvr (day-1)
  ! surface runoff (adjustable)
  REAL(WP)                             :: SAREAMAX    ! maximum saturated area
  REAL(WP)                             :: AXV_BEXP    ! ARNO/VIC "b" exponent
  REAL(WP)                             :: LOGLAMB     ! mean value of the log-transformed topographic index (m)
  REAL(WP)                             :: TISHAPE     ! shape parameter for the topo index Gamma distribution (-)
  ! time delay in runoff
  REAL(WP)                             :: TIMEDELAY   ! time delay in runoff (days)
  ! snow model
  REAL(WP)                             :: MBASE       ! base melt temperature (deg. C)
  REAL(WP)                             :: MFMAX       ! maximum melt factor (mm melt deg C.-1 6hrs-1)
  REAL(WP)                             :: MFMIN       ! minimum melt factor (mm melt deg C.-1 6hrs-1)
  REAL(WP)                             :: PXTEMP      ! rain-snow partition temperature (deg. C)
  REAL(WP)                             :: OPG         ! precipitation gradient (-) 
  REAL(WP)                             :: LAPSE       ! temperature gradient (deg. C)
 END TYPE PARADJ
 
 ! --------------------------------------------------------------------------------------
 ! (3) DERIVED PARAMETERS
 ! --------------------------------------------------------------------------------------
 TYPE PARDVD
  ! bucket sizes (derived)
  REAL(WP)                             :: MAXTENS_1   ! maximum tension storage in layer1 (mm)
  REAL(WP)                             :: MAXTENS_2   ! maximum tension storage in layer2 (mm)
  REAL(WP)                             :: MAXFREE_1   ! maximum free storage in layer 1 (mm)
  REAL(WP)                             :: MAXFREE_2   ! maximum free storage in layer2 (mm)
  REAL(WP)                             :: MAXTENS_1A  ! maximum storage in the recharge zone (mm)
  REAL(WP)                             :: MAXTENS_1B  ! maximum storage in the lower zone (mm)
  REAL(WP)                             :: MAXFREE_2A  ! maximum storage in the primary resvr (mm)
  REAL(WP)                             :: MAXFREE_2B  ! maximum storage in the secondary resvr (mm)
  ! evaporation
  REAL(WP)                             :: RTFRAC2     ! fraction of roots in the lower layer (-)
  ! percolation/baseflow
  REAL(WP)                             :: QBSAT       ! baseflow at saturation
  ! surface runoff
  REAL(WP)                             :: POWLAMB     ! mean value of the power-transformed topographic index (m**(1/n))
  REAL(WP)                             :: MAXPOW      ! max value of the power-transformed topographic index (m**(1/n))
  ! routing
  REAL(WP), DIMENSION(NTDH_MAX)        :: FRAC_FUTURE ! fraction of runoff in future time steps
  INTEGER(I4B)                         :: NTDH_NEED   ! number of time-steps with non-zero routing contribution
 END TYPE PARDVD
 
 ! --------------------------------------------------------------------------------------
 ! (4) LIST OF PARAMETERS FOR A GIVEN MODEL
 ! --------------------------------------------------------------------------------------
 TYPE PAR_ID
  CHARACTER(LEN=9)                     :: PARNAME     ! list of parameter names
 ENDTYPE PAR_ID
 
END MODULE multiparam_types
