module conv_funcs_module
USE nrtype                                 ! variable types
USE multiconst                             ! fixed parameters (lh vapzn, etc.)
implicit none
private
public::RELHM2SPHM,SPHM2RELHM,WETBULBTMP,rlhum2dewpt,dewpt2vpair
contains

! ----------------------------------------------------------------------
! series of functions to convert one thing to another
! (courtesy of Drew Slater)
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
FUNCTION MSLP2AIRP(MSLP, ELEV)
! compute air pressure using mean sea level pressure and elevation
! (after Shuttleworth, 1993)
!
! -- actually returns MSLP2AIRP in the same units as MSLP, because
!    ( (293.-0.0065*ELEV) / 293. )**5.256 is dimensionless
!
IMPLICIT NONE

REAL(WP), INTENT(IN)         :: MSLP      ! base pressure (Pa)
REAL(WP), INTENT(IN)         :: ELEV      ! elevation difference from base (m)

REAL(WP)                     :: MSLP2AIRP ! Air pressure (Pa)

MSLP2AIRP = MSLP * ( (293.-0.0065*ELEV) / 293. )**5.256

END FUNCTION MSLP2AIRP

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
FUNCTION RLHUM2DEWPT(T, RLHUM)
! Compute Dewpoint temperature from Relative Humidity
! ---- This is done with respect to water ONLY ----
!
! All units are SI standard - i.e. Kelvin and pascals
! Based on Tetens' formula (1930)
IMPLICIT NONE

REAL(WP), INTENT(IN)         :: T         ! Temperature           (K)
REAL(WP), INTENT(IN)         :: RLHUM     ! Relative Humidity     (%)


REAL(WP)                     :: RLHUM2DEWPT     ! Dewpoint Temp   (K)

REAL(WP)                     :: VPSAT     ! Sat. vapour pressure at T (Pa)
REAL(WP)                     :: TDCEL     ! Dewpoint temp Celcius (C)

! Units note :              Pa = N m-2 = kg m-1 s-2
! SATVPFRZ=     610.8       ! Saturation water vapour pressure at 273.16K (Pa)
! W_RATIO =       0.622     ! molecular weight ratio of water to dry air (-)

VPSAT = SATVPFRZ * EXP( (17.27*(T-TFREEZE)) / (237.30 + (T-TFREEZE)) ) ! sat vapor press at grid cell (Pa)
TDCEL = 237.30 * LOG( (VPSAT/SATVPFRZ)*(RLHUM/100.) ) / &              ! dewpoint temperature         (C)
        (17.27 - LOG( (VPSAT/SATVPFRZ)*(RLHUM/100.) ) ) 
RLHUM2DEWPT = TDCEL + TFREEZE

END FUNCTION RLHUM2DEWPT

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
FUNCTION DEWPT2RLHUM(T, DEWPT)
! Compute Relative humidity from dewpoint temperature
! ---- This is done with respect to water ONLY ----
!
! All units are SI standard - i.e. Kelvin and pascals
! Based on Tetens' formula (1930)
IMPLICIT NONE

REAL(WP), INTENT(IN)         :: T         ! Temperature           (K)
REAL(WP), INTENT(IN)         :: DEWPT     ! Dewpoint temp         (K)

REAL(WP)                     :: DEWPT2RLHUM ! Relative Humidity   (%)

REAL(WP)                     :: VPSAT     ! Sat. vapour pressure at T (Pa)
REAL(WP)                     :: TDCEL     ! Dewpt in celcius      (C)

! Units note :              Pa = N m-2 = kg m-1 s-2
! SATVPFRZ=     610.8       ! Saturation water vapour pressure at 273.16K (Pa)

TDCEL = DEWPT-TFREEZE
VPSAT = SATVPFRZ * EXP( (17.27*(T-TFREEZE)) / (237.30 + (T-TFREEZE)) )      ! Sat vapor press (Pa)
DEWPT2RLHUM = 100. * (SATVPFRZ/VPSAT) * EXP((17.27*TDCEL)/(237.30+TDCEL))   ! Relative Humidity (%)

END FUNCTION DEWPT2RLHUM


! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
FUNCTION DEWPT2SPHM(DEWPT, PRESS)
! Compute specific humidity from dewpoint temp with respect to water
! ---- This is done with respect to water ONLY ----
!
! All units are SI standard - i.e. Kelvin and pascals
! Based on Tetens' formula (1930)
! VPAIR is the current vapor pressure as it used dewpoint to compute staurated VP
IMPLICIT NONE

REAL(WP), INTENT(IN)         :: DEWPT     ! Dewpoint temp         (K)
REAL(WP), INTENT(IN)         :: PRESS     ! Pressure              (Pa)

REAL(WP)                     :: DEWPT2SPHM ! Specific Humidity    (g/g)

REAL(WP)                     :: VPAIR     ! vapour pressure at T  (Pa)
REAL(WP)                     :: TDCEL     ! Dewpt in celcius      (C)

! Units note :              Pa = N m-2 = kg m-1 s-2
! SATVPFRZ=     610.8       ! Saturation water vapour pressure at 273.16K (Pa)

TDCEL = DEWPT-TFREEZE
VPAIR = SATVPFRZ * EXP( (17.27*TDCEL) / (237.30 + TDCEL) )        ! Vapour Press           (Pa)
DEWPT2SPHM = (VPAIR * W_RATIO)/(PRESS - (1.-W_RATIO)*VPAIR)       ! Specific humidity (g/g)

END FUNCTION DEWPT2SPHM

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
FUNCTION DEWPT2VPAIR(DEWPT)
! Compute vapor pressure of the air from dewpoint temp with respect to water
! ---- This is done with respect to water ONLY ----
!
! All units are SI standard - i.e. Kelvin and pascals
! Based on Tetens' formula (1930)
! VPAIR is the current vapor pressure as it used dewpoint to compute staurated VP
IMPLICIT NONE

REAL(WP), INTENT(IN)         :: DEWPT     ! Dewpoint temp         (K)
REAL(WP)                     :: TDCEL     ! Dewpt in celcius      (C)

REAL(WP)                     :: DEWPT2VPAIR ! Vapour Press  (Pa)

! Units note :              Pa = N m-2 = kg m-1 s-2
! SATVPFRZ=     610.8       ! Saturation water vapour pressure at 273.16K (Pa)

TDCEL = DEWPT-TFREEZE
DEWPT2VPAIR = SATVPFRZ * EXP( (17.27*TDCEL) / (237.30 + TDCEL) )   ! Vapour Press  (Pa)

END FUNCTION DEWPT2VPAIR
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
FUNCTION SPHM2RELHM(SPHM, PRESS, TAIR)
! Compute specific humidity from dewpoint temp with respect to water
! ---- This is done with respect to water ONLY ----
!
! All units are SI standard - i.e. Kelvin and pascals
! Based on Tetens' formula (1930)
! VPAIR is the current vapor pressure as it used dewpoint to compute staurated VP
IMPLICIT NONE

REAL(WP), INTENT(IN)         :: SPHM      ! Specific Humidity (g/g)
REAL(WP), INTENT(IN)         :: PRESS     ! Pressure              (Pa)
REAL(WP), INTENT(IN)         :: TAIR      ! Air temp

REAL(WP)                     :: SPHM2RELHM ! Dewpoint Temp (K)

REAL(WP)                     :: VPSAT     ! vapour pressure at T  (Pa)
REAL(WP)                     :: TDCEL     ! Dewpt in celcius      (C)
!REAL(WP)                     :: DUM       ! Intermediate

! Units note :              Pa = N m-2 = kg m-1 s-2
! SATVPFRZ=     610.8       ! Saturation water vapour pressure at 273.16K (Pa)

TDCEL = TAIR-TFREEZE
VPSAT = SATVPFRZ * EXP( (17.27*TDCEL) / (237.30 + TDCEL) )       ! Vapour Press      (Pa)
SPHM2RELHM = (SPHM * PRESS)/(VPSAT * (W_RATIO + SPHM*(1.-W_RATIO)))

END FUNCTION SPHM2RELHM
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
FUNCTION RELHM2SPHM(RELHM, PRESS, TAIR)
! Compute specific humidity from dewpoint temp with respect to water
! ---- This is done with respect to water ONLY ----
!
! All units are SI standard - i.e. Kelvin and pascals
! Based on Tetens' formula (1930)
! VPAIR is the current vapor pressure as it used dewpoint to compute staurated VP
IMPLICIT NONE

REAL(WP), INTENT(IN)         :: RELHM     ! Relative Humidity     (%)
REAL(WP), INTENT(IN)         :: PRESS     ! Pressure              (Pa)
REAL(WP), INTENT(IN)         :: TAIR      ! Air temp

REAL(WP)                     :: RELHM2SPHM ! Specific Humidity (g/g)

REAL(WP)                     :: PVP       ! Partial vapour pressure at T  (Pa)
REAL(WP)                     :: TDCEL     ! Dewpt in celcius      (C)
!REAL(WP)                     :: DUM       ! Intermediate

! Units note :              Pa = N m-2 = kg m-1 s-2
! SATVPFRZ=     610.8       ! Saturation water vapour pressure at 273.16K (Pa)

TDCEL = TAIR-TFREEZE
PVP = RELHM * SATVPFRZ * EXP( (17.27*TDCEL)/(237.30 + TDCEL) ) ! Partial Vapour Press (Pa)
RELHM2SPHM = (PVP * W_RATIO)/(PRESS - (1. - W_RATIO)*PVP)

END FUNCTION RELHM2SPHM

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
FUNCTION WETBULBTMP(TAIR, RELHM, PRESS)
! Compute wet bulb temperature based on humidity and pressure
IMPLICIT NONE
! input
REAL(WP), INTENT(IN)         :: TAIR      ! Air temp              (K)
REAL(WP), INTENT(IN)         :: RELHM     ! Relative Humidity     (-)
REAL(WP), INTENT(IN)         :: PRESS     ! Pressure              (Pa)
! output
REAL(WP)                     :: WETBULBTMP ! Wet bulb temperature (K)
! locals
REAL(WP)                     :: Tcel           ! Temperature in celcius      (C)
REAL(WP)                     :: PVP            ! Partial vapor pressure (Pa)
REAL(WP)                     :: TWcel          ! Wet bulb temperature in celcius (C)
REAL(WP),PARAMETER           :: k=6.54E-4_WP   ! normalizing factor in wet bulb estimate (C-1)
REAL(WP)                     :: Twet_trial0    ! trial value for wet bulb temperature (C) 
REAL(WP)                     :: Twet_trial1    ! trial value for wet bulb temperature (C) 
REAL(WP)                     :: f0,f1          ! function evaluations (C)
REAL(WP)                     :: df_dT          ! derivative (-)
REAL(WP)                     :: TWinc          ! wet bulb temperature increment (C)
INTEGER(I4B)                 :: iter           ! iterattion index
REAL(WP),PARAMETER           :: Xoff=1.E-5_WP  ! finite difference increment (C)
REAL(WP),PARAMETER           :: Xtol=1.E-8_WP  ! convergence tolerance (C)
INTEGER(I4B)                 :: maxiter=15     ! maximum number of iterations 
! convert temperature to Celcius
Tcel = TAIR-TFREEZE
! compute partial vapor pressure based on temperature (Pa)
PVP = RELHM * SATVPRESS(Tcel)
! define an initial trial value for wetbulb temperature
TWcel = Tcel - 5._wp
! iterate until convergence
do iter=1,maxiter
 ! compute Twet estimates
 Twet_trial0 = Tcel - (SATVPRESS(TWcel)      - PVP)/(k*PRESS)
 Twet_trial1 = Tcel - (SATVPRESS(TWcel+Xoff) - PVP)/(k*PRESS)
 ! compute function evaluations
 f0 = Twet_trial0 - TWcel
 f1 = Twet_trial1 - (TWcel+Xoff)
 ! compute derivative and iteration increment
 df_dT = (f0 - f1)/Xoff
 TWinc = f0/df_dT
 ! compute new value of wet bulb temperature (C)
 TWcel = TWcel + TWinc
 ! check if achieved tolerance
 if(abs(f0) < Xtol) exit
 ! check convergence
 if(iter==maxiter)stop 'failed to converge in WETBULBTMP'
enddo  ! (iterating)

! return value in K
WETBULBTMP = TWcel + TFREEZE

END FUNCTION WETBULBTMP


! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
FUNCTION SATVPRESS(TCEL)
! Compute saturated vapor pressure (Pa)
! Units note :              Pa = N m-2 = kg m-1 s-2
! SATVPFRZ=     610.8       ! Saturation water vapour pressure at 273.16K (Pa)
IMPLICIT NONE
REAL(WP),INTENT(IN) :: TCEL      ! Temperature (C)
REAL(WP)            :: SATVPRESS ! Saturated vapor pressure (Pa)
SATVPRESS = SATVPFRZ * EXP( (17.27*TCEL)/(237.30 + TCEL) ) ! Saturated Vapour Press (Pa)
END FUNCTION SATVPRESS


end module conv_funcs_module
