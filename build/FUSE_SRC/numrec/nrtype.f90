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

  ! define string lengths
  integer(i4b), parameter :: strLen     = 256
  integer(i4b), parameter :: FileStrLen = 300
  integer(i4b), parameter :: gageStrLen = 30

END MODULE nrtype
