MODULE multiforce_types

 ! ---------------------------------------------------------------------------------------
 ! Creator:
 ! --------
 ! Martyn Clark
 ! Modified by Brian Henn to include snow model, 6/2013
 ! Modified by Nans Addor to enable distributed modeling, 9/2016
 ! Modified by Cyril Thébault to allow different metrics as objective function, 2024
 ! Modified by Martyn Clark to separate type definitions from data storage, 01/2026
 ! ---------------------------------------------------------------------------------------
 
 USE nrtype
 
 implicit none
 private

 public :: TDATA, VDATA, FDATA

 ! the time data structure (will have no spatial dimension)
 TYPE TDATA
    INTEGER(I4B)                         :: IY         ! year
    INTEGER(I4B)                         :: IM         ! month
    INTEGER(I4B)                         :: ID         ! day
    INTEGER(I4B)                         :: IH         ! hour
    INTEGER(I4B)                         :: IMIN       ! minute
    REAL(WP)                             :: DSEC       ! second
    REAL(WP)                             :: DTIME      ! time in seconds since year dot
 ENDTYPE TDATA

 ! validation/observation data
 !   gridded data:    (nx, ny, nt)
 !   point/HRU data:  (1, nobs, nt), with singleton first spatial dimension
 TYPE VDATA
    REAL(WP)                             :: OBSQ       ! observed runoff (mm day-1)
 END TYPE VDATA
 
 ! meteorological forcing data
 !   gridded data:    (nx, ny, nt)
 !   point/HRU data:  (1, nhru, nt), with singleton first spatial dimension
 TYPE FDATA
    REAL(WP)                             :: PPT        ! water input: rain + melt (mm day-1)
    REAL(WP)                             :: TEMP       ! temperature for snow model (deg.C)
    REAL(WP)                             :: PET        ! energy input: potential ET (mm day-1)
 ENDTYPE FDATA

END MODULE multiforce_types
