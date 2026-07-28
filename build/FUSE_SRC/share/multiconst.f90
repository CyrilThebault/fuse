MODULE multiconst
 USE nrtype
 ! define physical constants
 REAL(WP), PARAMETER         :: ave_slp      =  101325.0_wp      ! mean sea level pressure              (Pa)
 REAL(WP), PARAMETER         :: vkc          =       0.4_wp      ! von Karman constant                  (-)
 REAL(WP), PARAMETER         :: satvpfrz     =     610.8_wp      ! sat vapour pressure at 273.16K       (Pa)
 REAL(WP), PARAMETER         :: w_ratio      =       0.622_wp    ! molecular ratio water to dry air     (-)
 REAL(WP), PARAMETER         :: R_da         =     287.053_wp    ! gas constant for dry air             (Pa K-1 m3 kg-1; J kg-1 K-1)
 REAL(WP), PARAMETER         :: R_wv         =     461.285_wp    ! gas constant for water vapor         (Pa K-1 m3 kg-1; J kg-1 K-1)
 REAL(WP), PARAMETER         :: gravity      =       9.80616_wp  ! acceleration of gravity              (m s-2)
 REAL(WP), PARAMETER         :: Cp_air       =    1005._wp       ! specific heat of air                 (J kg-1 K-1)
 REAL(WP), PARAMETER         :: Cp_ice       =    2114._wp       ! specific heat of ice                 (J kg-1 K-1)
 REAL(WP), PARAMETER         :: Cp_soil      =     850._wp       ! specific heat of soil                (J kg-1 K-1)
 REAL(WP), PARAMETER         :: Cp_water     =    4181._wp       ! specific heat of liquid water        (J kg-1 K-1)
 REAL(WP), PARAMETER         :: Tfreeze      =     273.16_wp     ! temperature at freezing              (K)
 REAL(WP), PARAMETER         :: TriplPt      =     273.16_wp     ! triple point of water                (K)
 REAL(WP), PARAMETER         :: LH_fus       =  333700.0_wp      ! latent heat of fusion                (J kg-1)
 REAL(WP), PARAMETER         :: LH_vap       = 2501000.0_wp      ! latent heat of vaporization          (J kg-1)
 REAL(WP), PARAMETER         :: LH_sub       = 2834700.0_wp      ! latent heat of sublimation           (J kg-1)
 REAL(WP), PARAMETER         :: sigma        =       5.6705d-8   ! Stefan Boltzman constant             (W m-2 K-4)
 REAL(WP), PARAMETER         :: em_sno       =       0.99_wp     ! emissivity of snow                   (-)
 REAL(WP), PARAMETER         :: lambda_air   =       0.023_wp    ! thermal conductivity of air          (W m-1 K-1)
 REAL(WP), PARAMETER         :: lambda_ice   =       2.29_wp     ! thermal conductivity of ice          (W m-1 K-1)
 REAL(WP), PARAMETER         :: lambda_soil  =       3.21_wp     ! thermal conductivity of soil         (W m-1 K-1)
 REAL(WP), PARAMETER         :: lambda_water =       0.60_wp     ! thermal conductivity of liquid water (W m-1 K-1)
 REAL(WP), PARAMETER         :: iden_air     =       1.293_wp    ! intrinsic density of air             (kg m-3)
 REAL(WP), PARAMETER         :: iden_ice     =     917.0_wp      ! intrinsic density of ice             (kg m-3)
 REAL(WP), PARAMETER         :: iden_water   =    1000.0_wp      ! intrinsic density of liquid water    (kg m-3)
 REAL(WP), PARAMETER         :: secprday     =   86400._wp       ! number of seconds in a day
 REAL(WP), PARAMETER         :: secprhour    =    3600._wp       ! number of seconds in an hour
 REAL(WP), PARAMETER         :: secprmin     =      60._wp       ! number of seconds in a minute
END MODULE multiconst
