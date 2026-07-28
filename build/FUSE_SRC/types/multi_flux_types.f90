MODULE multi_flux_types

 USE nrtype

 implicit none
 private

 public :: FLUXES

 TYPE FLUXES
  REAL(WP)                             :: EFF_PPT     ! effective precipitation (mm day-1)
  REAL(WP)                             :: SATAREA     ! saturated area (-)
  REAL(WP)                             :: QSURF       ! surface runoff (mm day-1)
  REAL(WP)                             :: EVAP_1A     ! evaporation from soil excess zone (mm day-1)
  REAL(WP)                             :: EVAP_1B     ! evaporation from soil recharge zone (mm day-1)
  REAL(WP)                             :: EVAP_1      ! evaporation from upper soil layer (mm day-1)
  REAL(WP)                             :: EVAP_2      ! evaporation from lower soil layer (mm day-1)
  REAL(WP)                             :: RCHR2EXCS   ! flow from recharge to excess (mm day-1)
  REAL(WP)                             :: TENS2FREE_1 ! flow from tension storage to free storage (mm day-1)
  REAL(WP)                             :: TENS2FREE_2 ! flow from tension storage to free storage (mm day-1)
  REAL(WP)                             :: QINTF_1     ! interflow from free water (mm day-1)
  REAL(WP)                             :: QPERC_12    ! percolation from upper to lower soil layers (mm day-1)
  REAL(WP)                             :: QBASE_2     ! baseflow (mm day-1)
  REAL(WP)                             :: QBASE_2A    ! baseflow from primary linear resvr (mm day-1)
  REAL(WP)                             :: QBASE_2B    ! baseflow from secondary linear resvr (mm day-1)
  REAL(WP)                             :: OFLOW_1     ! bucket overflow (mm day-1)
  REAL(WP)                             :: OFLOW_2     ! bucket overflow (mm day-1)
  REAL(WP)                             :: OFLOW_2A    ! bucket overflow (mm day-1)
  REAL(WP)                             :: OFLOW_2B    ! bucket overflow (mm day-1)
  REAL(WP)                             :: ERR_WATR_1  ! excessive extrapolation: total storage in layer1 (mm day-1)
  REAL(WP)                             :: ERR_TENS_1  ! excessive extrapolation: tension storage in layer1 (mm day-1)
  REAL(WP)                             :: ERR_FREE_1  ! excessive extrapolation: free storage in layer 1 (mm day-1)
  REAL(WP)                             :: ERR_TENS_1A ! excessive extrapolation: storage in the recharge zone (mm day-1)
  REAL(WP)                             :: ERR_TENS_1B ! excessive extrapolation: storage in the lower zone (mm day-1)
  REAL(WP)                             :: ERR_WATR_2  ! excessive extrapolation: total storage in layer2 (mm day-1)
  REAL(WP)                             :: ERR_TENS_2  ! excessive extrapolation: tension storage in layer2 (mm day-1)
  REAL(WP)                             :: ERR_FREE_2  ! excessive extrapolation: free storage in layer2 (mm day-1)
  REAL(WP)                             :: ERR_FREE_2A ! excessive extrapolation: storage in the primary resvr (mm day-1)
  REAL(WP)                             :: ERR_FREE_2B ! excessive extrapolation: storage in the secondary resvr (mm day-1)
  REAL(WP)                             :: CHK_TIME    ! time elapsed during time step (days)
 ENDTYPE FLUXES

END MODULE multi_flux_types
