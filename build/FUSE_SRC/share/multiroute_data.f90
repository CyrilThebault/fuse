MODULE multiroute

 USE nrtype
 USE multiroute_types, only: RUNOFF

 implicit none
 private

 public :: MROUTE

 TYPE(RUNOFF)                                 :: MROUTE     ! runoff for one time step

END MODULE multiroute
