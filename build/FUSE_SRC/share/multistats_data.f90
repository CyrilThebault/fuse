MODULE multistats

 USE nrtype
 USE multistats_types, only: SUMMARY

 implicit none
 private

 public :: MSTATS, MOD_IX, FCOUNT
 
 TYPE(SUMMARY)                         :: MSTATS        ! (model summary statistics)
 INTEGER(I4B)                          :: MOD_IX = 1    ! (model index)
 INTEGER(I4B)                          :: FCOUNT        ! (number of model simulations)

END MODULE multistats
