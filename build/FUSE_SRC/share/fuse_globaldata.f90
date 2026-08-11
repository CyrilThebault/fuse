MODULE fuse_globaldata

 USE nrtype
 
 implicit none
 include "fuseversion.inc"

 ! time step
 REAL(WP), save          :: CURRENT_DT            ! current time step (days)

 ! missing values
 INTEGER(I4B),PARAMETER  :: NA_VALUE=-9999        ! integer designating missing values - TODO: retrieve from NetCDF file
 REAL(WP),PARAMETER      :: NA_VALUE_SP=-9999_wp  ! integer designating missing values - TODO: retrieve from NetCDF file

 ! NetCDF
 integer(i4b), save      :: ncid_out=-1           ! NetCDF output file ID

 ! initial store fraction (initialization)
 real(wp), parameter     :: fracState0=0.25_wp

 ! original code
 logical(lgt), save      :: isOriginal=.true.

 ! print flag
 logical(lgt), save      :: isPrint=.true.
 logical(lgt), save      :: isDebug=.false.

 ! flag to run mizuRoute
 logical(lgt), save      :: do_mizuRoute=.false.

 ! indices of forcing variables in vectors
 integer(i4b), parameter :: NVAR_HYDROMET=4
 integer(i4b), parameter :: iPRECIP=1, iTEMP=2, iPET=3, iQOBS=4

 ! indices of snow parameters in vectors
 integer(i4b), parameter :: NPAR_SNOW=7
 integer(i4b), parameter :: iMBASE=1, iMFMAX=2, iMFMIN=3, iPXTEMP=4, iOPG=5, iLAPSE=6  ! indices in vectors
 integer(i4b), parameter :: iPERR=7   ! not a snow parameter, but used here

 ! number of fuse evaluations
 integer(i4b), save      :: nFUSE_eval

end MODULE fuse_globaldata
