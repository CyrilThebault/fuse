MODULE multistate_types

 USE nrtype

 implicit none
 private

 public :: STATEV, M_TIME

 ! --------------------------------------------------------------------------------------
 ! model state structure
 ! --------------------------------------------------------------------------------------
 TYPE STATEV
  ! interception layer
  REAL(WP)                             :: SINT_0     ! interception storage (mm)
  ! snow layer
  REAL(WP)                             :: SWE_TOT    ! total storage as snow (mm)
  ! upper layer
  REAL(WP)                             :: WATR_1     ! total storage in layer1 (mm)
  REAL(WP)                             :: TENS_1     ! tension storage in layer1 (mm)
  REAL(WP)                             :: FREE_1     ! free storage in layer 1 (mm)
  REAL(WP)                             :: TENS_1A    ! storage in the recharge zone (mm)
  REAL(WP)                             :: TENS_1B    ! storage in the lower zone (mm)
  ! lower layer
  REAL(WP)                             :: WATR_2     ! total storage in layer2 (mm)
  REAL(WP)                             :: TENS_2     ! tension storage in layer2 (mm)
  REAL(WP)                             :: FREE_2     ! free storage in layer2 (mm)
  REAL(WP)                             :: FREE_2A    ! storage in the primary resvr (mm)
  REAL(WP)                             :: FREE_2B    ! storage in the secondary resvr (mm)
 END TYPE STATEV
 
 ! --------------------------------------------------------------------------------------
 ! model time structure
 ! --------------------------------------------------------------------------------------
 TYPE M_TIME
  REAL(WP)                             :: STEP       ! (time interval to advance model states)
 END TYPE M_TIME

END MODULE multistate_types
