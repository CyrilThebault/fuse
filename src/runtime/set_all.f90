module set_all_module
USE nrtype
USE netcdf
implicit none

contains

  SUBROUTINE SET_STATE(VAL, FSTATE, MBANDS)

    ! ---------------------------------------------------------------------------------------
    ! Creator:
    ! --------
    ! Nans Addor based on Martyn Clark's INIT_STATE
    ! ---------------------------------------------------------------------------------------
    ! Purpose:
    ! --------
    ! Set model states to a given value - useful to set them to _FillValue
    ! ---------------------------------------------------------------------------------------

    use multistate_types,  only: STATEV
    use multibands_types,  only: BANDS_VAR

    IMPLICIT NONE
    REAL(WP)       , INTENT(IN)           :: VAL          ! value
    type(STATEV)   , intent(inout)        :: FSTATE       ! model states
    type(BANDS_VAR), intent(inout)        :: MBANDS(:)    ! elevation bands
    
    ! interception
    FSTATE%SINT_0  = VAL
    
    ! upper layer
    FSTATE%TENS_1A = VAL
    FSTATE%TENS_1B = VAL
    FSTATE%TENS_1  = VAL
    FSTATE%FREE_1  = VAL
    FSTATE%WATR_1  = VAL
    
    ! lower layer
    FSTATE%TENS_2  = VAL
    FSTATE%FREE_2  = VAL
    FSTATE%FREE_2A = VAL
    FSTATE%FREE_2B = VAL
    FSTATE%WATR_2  = VAL
    
    ! snow model
    MBANDS(:)%SWE  = VAL
    FSTATE%SWE_TOT = VAL

    ! ---------------------------------------------------------------------------------------
  END SUBROUTINE SET_STATE

  SUBROUTINE SET_FLUXES(VAL, M_FLUX, MBANDS)

    ! ---------------------------------------------------------------------------------------
    ! Creator:
    ! --------
    ! Nans Addor based on Martyn Clark's INITFLUXES
    ! ---------------------------------------------------------------------------------------
    ! Purpose:
    ! --------
    ! Set model fluxes to a given value - useful to set them to _FillValue
    ! ---------------------------------------------------------------------------------------

    use multi_flux_types,  only: FLUXES
    use multibands_types,  only: BANDS_VAR

    IMPLICIT NONE
    REAL(WP)       , INTENT(IN)           :: VAL          ! value
    type(FLUXES)   , intent(inout)        :: M_FLUX       ! model fluxes
    type(BANDS_VAR), intent(inout)        :: MBANDS(:)    ! elevation bands
    
    M_FLUX%EVAP_0      = VAL
    M_FLUX%PIN0        = VAL
    M_FLUX%PTHRU       = VAL
    M_FLUX%EFF_PPT     = VAL
    M_FLUX%SATAREA     = VAL
    M_FLUX%QSURF       = VAL
    M_FLUX%EVAP_1A     = VAL
    M_FLUX%EVAP_1B     = VAL
    M_FLUX%EVAP_1      = VAL
    M_FLUX%EVAP_2      = VAL
    M_FLUX%RCHR2EXCS   = VAL
    M_FLUX%TENS2FREE_1 = VAL
    M_FLUX%TENS2FREE_2 = VAL
    M_FLUX%QINTF_1     = VAL
    M_FLUX%QPERC_12    = VAL
    M_FLUX%QBASE_2     = VAL
    M_FLUX%QBASE_2A    = VAL
    M_FLUX%QBASE_2B    = VAL
    M_FLUX%OFLOW_1     = VAL
    M_FLUX%OFLOW_2     = VAL
    M_FLUX%OFLOW_2A    = VAL
    M_FLUX%OFLOW_2B    = VAL
    
    MBANDS(:)%SNOWACCMLTN = VAL
    MBANDS(:)%SNOWMELT    = VAL
    
    M_FLUX%ERR_WATR_1  = VAL
    M_FLUX%ERR_TENS_1  = VAL
    M_FLUX%ERR_FREE_1  = VAL
    M_FLUX%ERR_TENS_1A = VAL
    M_FLUX%ERR_TENS_1B = VAL
    M_FLUX%ERR_WATR_2  = VAL
    M_FLUX%ERR_TENS_2  = VAL
    M_FLUX%ERR_FREE_2  = VAL
    M_FLUX%ERR_FREE_2A = VAL
    M_FLUX%ERR_FREE_2B = VAL
    M_FLUX%CHK_TIME    = VAL

  END SUBROUTINE SET_FLUXES

  SUBROUTINE SET_ROUTE(VAL, MROUTE)

    ! ---------------------------------------------------------------------------------------
    ! Creator:
    ! --------
    ! Nans Addor based on Martyn Clark's INIT_STATE
    ! ---------------------------------------------------------------------------------------
    ! Purpose:
    ! --------
    ! Set runoff variables to a given value - useful to set them to _FillValue
    ! ---------------------------------------------------------------------------------------

    use multiroute_types,  only: RUNOFF

    IMPLICIT NONE
    REAL(WP)       , INTENT(IN)           :: VAL          ! value
    type(RUNOFF)   , intent(inout)        :: MROUTE       ! model routing structure
    
    MROUTE%Q_INSTNT = VAL     ! instantaneous runoff
    MROUTE%Q_ROUTED = VAL     ! routed runoff
    MROUTE%Q_ACCURATE  = VAL  ! "accurate" runoff estimate (mm day-1)

  END SUBROUTINE SET_ROUTE

  SUBROUTINE SET_SNOW(VAL, MBANDS)

    ! ---------------------------------------------------------------------------------------
    ! Creator:
    ! --------
    ! Nans Addor based on Martyn Clark's INIT_STATE
    ! ---------------------------------------------------------------------------------------
    ! Purpose:
    ! --------
    ! Set snow variables to a given value - useful to set them to _FillValue
    ! ---------------------------------------------------------------------------------------
    
    use multibands_types,  only: BANDS_VAR

    IMPLICIT NONE
    REAL(WP)       , INTENT(IN)           :: VAL          ! value
    type(BANDS_VAR), intent(inout)        :: MBANDS(:)    ! elevation bands

    MBANDS(:)%SWE         = VAL       ! band snowpack water equivalent (mm)
    MBANDS(:)%SNOWACCMLTN = VAL       ! new snow accumulation in band (mm day-1)
    MBANDS(:)%SNOWMELT    = VAL       ! snowmelt in band (mm day-1)
    MBANDS(:)%DSWE_DT     = VAL       ! rate of change of band SWE (mm day-1)

  END SUBROUTINE SET_SNOW

end module set_all_module
