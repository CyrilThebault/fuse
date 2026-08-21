MODULE multibands
 
 ! Created by Brian Henn to allow multi-band snow modeling, 6/2013
 ! Based on module MULTIFORCE by Martyn Clark
 
 ! Modified by Martyn Clark to separate type definitions from data storage, 01/2026

 USE nrtype

 USE multibands_types, only: BANDS

 implicit none
 private

 public :: N_BANDS
 public :: MBANDS
 public :: Z_FORCING

 ! --------------------------------------------------------------------------------------
 TYPE(BANDS),DIMENSION(:),ALLOCATABLE  :: MBANDS          ! basin band information

 INTEGER(I4B)                          :: N_BANDS=0       ! number of bands, initialize to zero
 REAL(WP)                              :: Z_FORCING       ! elevation of forcing data (m)
 ! --------------------------------------------------------------------------------------

END MODULE multibands
