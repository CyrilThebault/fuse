SUBROUTINE UPDATE_SWE(DT)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Brian Henn, as part of FUSE snow model implementation, 6/2013
! Based on subroutines QSATEXCESS and UPDATSTATE, by Martyn Clark
! Modified by Nans Addor to enable distributed modeling, 9/2016
! Modified by Martyn Clark to enable the split info/var structure, 01/2026
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Computes the snow accumulation and melt from forcing data
! Then updates the SWE band states based on the fluxes
! ---------------------------------------------------------------------------------------
! Modules Modified:
! -----------------
! MODULE multibands -- SWE bands stored in MODULE multibands
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE model_defn                                        ! model definition structure
USE model_defnames                                    ! integer model definitions
USE multiparam                                        ! model parameters
USE multibands                                        ! model basin band structure
USE multiforce                                        ! model forcing
USE multistate                                        ! model states
USE multi_flux                                        ! model fluxes
IMPLICIT NONE
! input
REAL(SP), INTENT(IN)                   :: DT             ! length of the time step
! internal variables
LOGICAL(LGT)                           :: LEAP           ! leap year flag
REAL(SP)                               :: JDAY           ! Julian day of year
REAL(SP), DIMENSION(12)                :: CUMD           ! cumulative days of year
REAL(SP)                               :: DZ             ! vert. distance from forcing
REAL(SP)                               :: TEMP_Z         ! band temperature at timestep
REAL(SP)                               :: PRECIP_Z       ! band precipitation at timestep
REAL(SP)                               :: MF             ! melt factor (mm/deg.C-6hr)
INTEGER(I4B)                           :: ISNW           ! loop through snow model bands
! ---------------------------------------------------------------------------------------
! snow accumulation and melt calculations for each band
! also calculates effective precipitation
! ---------------------------------------------------------------------------------------
! first calculate day of year for melt factor calculation
CUMD = real((/ 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 /),sp)
IF (MOD(timDat%IY,4).EQ.0) THEN
 LEAP = .TRUE.
 CUMD(3:12) = CUMD(3:12) + 1._sp
ELSE
 LEAP = .FALSE.
ENDIF

JDAY = CUMD(timDat%IM) + timDat%ID

IF (LEAP) THEN   ! calculate melt factor from MFMAX, MFMIN and day of year
 MF = ((0.5_sp*SIN(((JDAY-81._sp)*2._sp*PI)/366._sp))+0.5_sp)*(MPARAM%MFMAX - MPARAM%MFMIN) + MPARAM%MFMIN
ELSE
 MF = ((0.5_sp*SIN(((JDAY-80._sp)*2._sp*PI)/365._sp))+0.5_sp)*(MPARAM%MFMAX - MPARAM%MFMIN) + MPARAM%MFMIN
ENDIF

! loop through model bands
DO ISNW=1,N_BANDS

 ! ---------------------------------------------------------------------------------------
 associate( &    ! link to the info and var sub-structures in MBANDS (less invasive / more readable in code below)
    z_mid       => mbands(isnw)%info%z_mid,      &
    af          => mbands(isnw)%info%af,         &
    swe         => mbands(isnw)%var%swe,         &
    snowaccmltn => mbands(isnw)%var%snowaccmltn, &
    snowmelt    => mbands(isnw)%var%snowmelt,    &
    dswe_dt     => mbands(isnw)%var%dswe_dt )
 
  ! calculate forcing data for each band
 DZ = Z_MID - Z_FORCING
 TEMP_Z = MFORCE%TEMP + DZ*MPARAM%LAPSE/1000._sp ! adjust for elevation using lapse rate
 IF (DZ.GT.0._sp) THEN ! adjust for elevation using OPG
  PRECIP_Z = MFORCE%PPT * (1._sp + DZ*MPARAM%OPG/1000._sp)
 ELSE
  PRECIP_Z = MFORCE%PPT / (1._sp - DZ*MPARAM%OPG/1000._sp)
 ENDIF
 IF ((SWE.GT.0._sp).AND.(TEMP_Z.GT.MPARAM%MBASE)) THEN
  ! calculate the initial snowmelt rate from the melt factor and the temperature
  SNOWMELT = MF*(TEMP_Z - MPARAM%MBASE) ! MBANDS%SNOWMELT has units of mm day-1
 ELSE
  SNOWMELT = 0.0_sp
 ENDIF

 ! calculate the accumulation rate from the forcing data
 IF (TEMP_Z.LT.MPARAM%PXTEMP) THEN
  SELECT CASE(SMODL%iRFERR)
   CASE(iopt_additive_e) ! additive rainfall error
    SNOWACCMLTN = MAX(0.0_sp, PRECIP_Z + MPARAM%RFERR_ADD)
   CASE(iopt_multiplc_e) ! multiplicative rainfall error
    SNOWACCMLTN = PRECIP_Z * MPARAM%RFERR_MLT
   CASE DEFAULT       ! check for errors
    print *, "SMODL%iRFERR must be either iopt_additive_e or iopt_multiplc_e"
    STOP
  END SELECT
 ELSE
  SNOWACCMLTN = 0.0_sp
 ENDIF

 ! update SWE, and check to ensure non-negative values
 DSWE_DT = SNOWACCMLTN - SNOWMELT
 IF ((SWE + DSWE_DT*DT).GE.0._sp) THEN
  SWE = SWE + DSWE_DT*DT
 ELSE ! reduce melt rate in case of negative SWE
  SNOWMELT = SWE/DT + SNOWACCMLTN
  SWE = 0.0_sp
 ENDIF

 ! calculate rainfall plus snowmelt
 IF (TEMP_Z.GT.MPARAM%PXTEMP) THEN
  SELECT CASE(SMODL%iRFERR)
   CASE(iopt_additive_e) ! additive rainfall error
   M_FLUX%EFF_PPT = M_FLUX%EFF_PPT + AF * &
   (MAX(0.0_sp, PRECIP_Z + MPARAM%RFERR_ADD) + SNOWMELT)
   CASE(iopt_multiplc_e) ! multiplicative rainfall error
   M_FLUX%EFF_PPT = M_FLUX%EFF_PPT + AF * &
   (PRECIP_Z * MPARAM%RFERR_MLT +  SNOWMELT)
   CASE DEFAULT       ! check for errors
    print *, "SMODL%iRFERR must be either iopt_additive_e or iopt_multiplc_e"
    STOP
  END SELECT
 ELSE
  M_FLUX%EFF_PPT = M_FLUX%EFF_PPT + AF * SNOWMELT
 ENDIF

 end associate

END DO  ! looping through bands

END SUBROUTINE UPDATE_SWE
