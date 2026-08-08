MODULE multiroute

 USE nrtype
 USE multiroute_types, only: RUNOFF

 implicit none
 private

 public :: FUTURE
 public :: MROUTE

 REAL(WP), ALLOCATABLE                        :: FUTURE(:)  ! runoff placed in future time steps

 TYPE(RUNOFF)                                 :: MROUTE     ! runoff for one time step

END MODULE multiroute
