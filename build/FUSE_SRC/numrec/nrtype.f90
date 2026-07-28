MODULE nrtype
  INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
  INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(4)
  INTEGER, PARAMETER :: I1B = SELECTED_INT_KIND(2)
  INTEGER, PARAMETER :: DP  = KIND(1.d0)
  INTEGER, PARAMETER :: SP  = KIND(1.0)
  INTEGER, PARAMETER :: SPC = KIND((1.0,1.0))
  INTEGER, PARAMETER :: LGT = KIND(.true.)
  integer, parameter :: I8B = selected_int_kind(15)
  integer, parameter :: DPC = kind((1.0D0,1.0D0))

  ! define specific precision
  integer, parameter :: WP  = DP ! working precision
  INTEGER, PARAMETER :: WPC = KIND((1.0_wp,1.0_wp))
  INTEGER, PARAMETER :: MSP = KIND(1.0)  ! SP still needed for f77 netcdf routines

  ! define constants
  REAL(WP), PARAMETER :: PI_SP=3.141592653589793238462643383279502884197_wp
  REAL(WP), PARAMETER :: PIO2=1.57079632679489661923132169163975144209858_wp
  REAL(WP), PARAMETER :: TWOPI=6.283185307179586476925286766559005768394_wp
  REAL(WP), PARAMETER :: SQRT2=1.41421356237309504880168872420969807856967_wp
  REAL(WP), PARAMETER :: EULER=0.5772156649015328606065120900824024310422_wp

  ! define string lengths
  integer(i4b), parameter :: strLen     = 256
  integer(i4b), parameter :: FileStrLen = 300
  integer(i4b), parameter :: gageStrLen = 30

  TYPE sprs2_wp
    INTEGER(I4B) :: n,len
    REAL(WP), DIMENSION(:), POINTER :: val
    INTEGER(I4B), DIMENSION(:), POINTER :: irow
    INTEGER(I4B), DIMENSION(:), POINTER :: jcol
  END TYPE sprs2_wp
!   TYPE sprs2_wp
!     INTEGER(I4B) :: n,len
!     REAL(WP), DIMENSION(:), POINTER :: val
!     INTEGER(I4B), DIMENSION(:), POINTER :: irow
!     INTEGER(I4B), DIMENSION(:), POINTER :: jcol
!   END TYPE sprs2_wp
END MODULE nrtype
