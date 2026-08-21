SUBROUTINE MEAN_TIPOW()
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2007
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Computes the mean of the power-transformed topographic index
! ---------------------------------------------------------------------------------------
! Modules Modified:
! -----------------
! MODULE multiparam -- mean topographic index stored in MODULE multiparam
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE nr, ONLY : gammp                                  ! interface for the incomplete gamma function
USE multiparam                                        ! model parameters
IMPLICIT NONE
! internal variables
INTEGER(I4B)                           :: IBIN        ! loop through bins
INTEGER(I4B), PARAMETER                :: NBINS=2000  ! number of bins in PDF of topo index
REAL(WP), PARAMETER                    :: TI_MAX=50._WP  ! maximum possible log-transformed index
REAL(WP)                               :: TI_OFF      ! offset in the Gamma distribution
REAL(WP)                               :: TI_SHP      ! shape of the Gamma distribution
REAL(WP)                               :: TI_CHI      ! CHI, see Sivapalan et al., 1987
REAL(WP)                               :: LOWERV      ! lower value of frequency bin
REAL(WP)                               :: UPPERV      ! upper value of frequency bin
REAL(WP)                               :: LOWERP      ! cumulative probability of the lower value
REAL(WP)                               :: UPPERP      ! cumulative probability of the upper value
REAL(WP)                               :: GMARG2      ! 2nd argument to the incomplete Gamma function
REAL(WP)                               :: PROBIN      ! probability of the current bin
REAL(WP)                               :: LOGVAL      ! log-transformed index for the current bin
REAL(WP)                               :: POWVAL      ! power-transformed index for the current bin
!REAL(WP)                               :: AVELOG      ! average log-transformed index (testing)
REAL(WP)                               :: AVEPOW      ! average power-transformed index
! ---------------------------------------------------------------------------------------
! preliminaries -- get parameters of the Gamma distribution (save typing)
TI_OFF = 3._WP              ! offset in the Gamma distribution (the "3rd" parameter)
TI_SHP = MPARAM%TISHAPE  ! shape of the Gamma distribution (the "2nd" parameter)
TI_CHI = (MPARAM%LOGLAMB - TI_OFF) / MPARAM%TISHAPE ! Chi -- loglamb is the first parameter (mean)
! values for testing (Sivapalan et al., WRR, December 1987)
!TI_OFF = 3.82_WP  ! TI_OFF = 2.92_WP
!TI_SHP = 2.48_WP  ! TI_SHP = 3.52_WP
!TI_CHI = 1.00_WP  ! TI_CHI = 0.742_WP
! loop through the frequency distribution
LOWERV = 0._WP
LOWERP = 0._WP
!AVELOG = 0._WP
AVEPOW = 0._WP
DO IBIN=1,NBINS
 ! get probability for the current bin
 UPPERV = (REAL(IBIN)/REAL(NBINS)) * TI_MAX          ! upper value in frequency bin
 GMARG2 = MAX(0._WP, UPPERV - TI_OFF) / TI_CHI          ! 2nd argument to the Gamma function
 UPPERP = GAMMP(TI_SHP, GMARG2)                      ! GAMMP is the incomplete Gamma function
 PROBIN = UPPERP-LOWERP                              ! probability of the current bin
 ! get the scaled topographic index value
 LOGVAL = 0.5_WP*(LOWERV+UPPERV)                        ! log-transformed index for the current bin
 POWVAL = (EXP(LOGVAL))**(1._WP/MPARAM%QB_POWR)         ! power-transformed index for the current bin
 !AVELOG = AVELOG + LOGVAL*PROBIN                     ! average log-transformed index (testing)
 AVEPOW = AVEPOW + POWVAL*PROBIN                     ! average power-transformed index
 !write(*,'(7(f9.3,1x))') lowerv, upperv, logval, powval, avelog, avepow
 ! save the lower value and probability
 LOWERV = UPPERV                                     ! lower value for the next bin
 LOWERP = UPPERP                                     ! cumulative probability for the next bin
END DO  ! (looping through bins)
DPARAM%MAXPOW  = POWVAL
DPARAM%POWLAMB = AVEPOW
!print *, DPARAM%POWLAMB, MPARAM%QB_POWR
!pause
! ---------------------------------------------------------------------------------------
END SUBROUTINE MEAN_TIPOW
