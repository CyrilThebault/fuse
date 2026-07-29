module sce_driver_MODULE

  USE nrtype
  use info_types,  only: fuse_info
  use work_types,  only: fuse_work
  use data_types,  only: domain_data

  use sce_callback_context, only: set_sce_context, clear_sce_context

  implicit none

  private
  public :: sce_driver

contains

  subroutine sce_driver(info, work, domain, APAR, BL, BU)
  USE multiparam, only: MAXN    ! maximum number of trials before optimization is terminated
  USE multiparam, only: KSTOP   ! number of shuffling loops the value must change by PCENTO
  USE multiparam, only: PCENTO  ! the percentage
  USE multiparam, only: NUMPAR  ! # parameters
  USE multiparam, only: LPARAM, PARATT

  USE GETPAR_STR_MODULE, only: GETPAR_STR

  USE parameter_transform_module, only: validate_transform
  USE parameter_transform_module, only: vector_to_search_space

  USE fuse_globaldata, only: isPrint ! used to turn of printing for calibration runs
  USE fuse_globaldata, only: nFUSE_eval ! # FUSE evaluations
  USE model_defn, only: FNAME_TEMPRY, FNAME_ASCII
  implicit none
  ! input/output
  type(fuse_info)       , intent(inout)  :: info     ! info structures (runtime settings etc.)
  type(fuse_work)       , intent(inout)  :: work     ! work structures that depend on npar/nState
  type(domain_data)     , intent(inout)  :: domain   ! the fuse domain structure that stores data arrays
  real(wp)              , intent(in)     :: APAR(:)  ! model parameter set
  real(wp)              , intent(in)     :: BL(:)    ! vector of lower parameter bounds
  real(wp)              , intent(in)     :: BU(:)    ! vector of upper parameter bounds
  ! internal variables
  REAL(MSP)                              :: AF_MSP    ! objective function value
  REAL(MSP), DIMENSION(:), ALLOCATABLE   :: APAR_MSP  ! ! lower bound of model parameters
  REAL(MSP), DIMENSION(:), ALLOCATABLE   :: BL_MSP    ! ! lower bound of model parameters
  REAL(MSP), DIMENSION(:), ALLOCATABLE   :: BU_MSP    ! ! upper bound of model parameters

  REAL(MSP), DIMENSION(:), ALLOCATABLE   :: APAR_PHYS_MSP
  REAL(MSP), DIMENSION(:), ALLOCATABLE   :: BL_PHYS_MSP
  REAL(MSP), DIMENSION(:), ALLOCATABLE   :: BU_PHYS_MSP

  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: TRANSFORM_CODES

  TYPE(PARATT)                           :: PARAM_META
  INTEGER(I4B)                           :: IPAR
  INTEGER(I4B)                           :: IERR
  CHARACTER(LEN=256)                     :: MESSAGE

  INTEGER(I4B)                           :: NOPT    ! number of parameters to be optimized
  INTEGER(I4B)                           :: NGS     ! # complexes in the initial population
  INTEGER(I4B)                           :: NPG     ! # points in each complex
  INTEGER(I4B)                           :: NPS     ! # points in a sub-complex
  INTEGER(I4B)                           :: NSPL    ! # evolution steps allowed for each complex before shuffling
  INTEGER(I4B)                           :: MINGS   ! minimum number of complexes required
  INTEGER(I4B)                           :: INIFLG  ! 1 = include initial point in the population
  INTEGER(I4B)                           :: IPRINT  ! 0 = supress printing
  INTEGER(I4B)                           :: ISCE    ! unit number for SCE write
  INTEGER(KIND=4)                        :: ISEED   ! seed for the random sequence

  NOPT   =  NUMPAR         ! number of parameters to be optimized (NUMPAR in module multiparam)
  NGS    =     10          ! number of complexes in the initial population
  NPG    =  2*NOPT + 1     ! number of points in each complex
  NPS    =    NOPT + 1     ! number of points in a sub-complex
  NSPL   =  2*NOPT + 1     ! number of evolution steps allowed for each complex before shuffling
  MINGS  =  NGS            ! minimum number of complexes required
  INIFLG =  1              ! 1 = include initial point in the population
  IPRINT =  1              ! 0 = supress printing

  ! Store physical-space values using the precision required by SCE.
  ALLOCATE(APAR_PHYS_MSP(NUMPAR))
  ALLOCATE(BL_PHYS_MSP(NUMPAR))
  ALLOCATE(BU_PHYS_MSP(NUMPAR))

  ALLOCATE(APAR_MSP(NUMPAR))
  ALLOCATE(BL_MSP(NUMPAR))
  ALLOCATE(BU_MSP(NUMPAR))

  ALLOCATE(TRANSFORM_CODES(NUMPAR))

  APAR_PHYS_MSP = APAR
  BL_PHYS_MSP   = BL
  BU_PHYS_MSP   = BU

  ! Retrieve and validate the transformation code for each model parameter.
  DO IPAR = 1, NUMPAR

    CALL GETPAR_STR(LPARAM(IPAR)%PARNAME, PARAM_META)

    TRANSFORM_CODES(IPAR) = PARAM_META%PARVTN

    CALL validate_transform(                  &
      LPARAM(IPAR)%PARNAME,                   &
      BL_PHYS_MSP(IPAR),                      &
      APAR_PHYS_MSP(IPAR),                    &
      BU_PHYS_MSP(IPAR),                      &
      TRANSFORM_CODES(IPAR),                  &
      IERR,                                   &
      MESSAGE)

    IF (IERR /= 0) THEN
      WRITE(*,'(A)') TRIM(MESSAGE)
      STOP 'Invalid parameter transformation'
    END IF

  END DO

  ! Transform the initial parameter set into optimizer search space.
  CALL vector_to_search_space(                &
    APAR_PHYS_MSP,                            &
    TRANSFORM_CODES,                          &
    APAR_MSP,                                 &
    IERR,                                     &
    MESSAGE)

  IF (IERR /= 0) THEN
    WRITE(*,'(A)') TRIM(MESSAGE)
    STOP 'Unable to transform initial parameter set'
  END IF

  ! Transform the lower parameter bounds into optimizer search space.
  CALL vector_to_search_space(                &
    BL_PHYS_MSP,                              &
    TRANSFORM_CODES,                          &
    BL_MSP,                                   &
    IERR,                                     &
    MESSAGE)

  IF (IERR /= 0) THEN
    WRITE(*,'(A)') TRIM(MESSAGE)
    STOP 'Unable to transform lower parameter bounds'
  END IF

  ! Transform the upper parameter bounds into optimizer search space.
  CALL vector_to_search_space(                &
    BU_PHYS_MSP,                              &
    TRANSFORM_CODES,                          &
    BU_MSP,                                   &
    IERR,                                     &
    MESSAGE)

  IF (IERR /= 0) THEN
    WRITE(*,'(A)') TRIM(MESSAGE)
    STOP 'Unable to transform upper parameter bounds'
  END IF

  ! pass the FUSE structures to the context setter
  ! NOTE: in sce_context_set, info/work/domain have the target attribute so can point to them
  call set_sce_context(info, work, domain, TRANSFORM_CODES)

  ! open up ASCII output file
  ISCE = 96 ! (file unit)
  FNAME_ASCII = trim(FNAME_TEMPRY)//'_sce_output.txt'
  print *, 'Creating SCE output file:', trim(FNAME_ASCII)
  OPEN(96, FILE=TRIM(FNAME_ASCII) )

  ! printing
  isPrint     = .false.  ! turn off printing to screen
  nFUSE_eval  = 0        ! number of fuse evaluations

  ! set random seed
  ISEED = 1

  ! optimize (returns A and AF)
  ! note that SCE requires the kind of APAR, BL, BU to be MSP
  CALL SCEUA(APAR_MSP,AF_MSP,BL_MSP,BU_MSP,NOPT,MAXN,KSTOP,PCENTO,ISEED,&
             NGS,NPG,NPS,NSPL,MINGS,INIFLG,IPRINT,ISCE)

  ! close ASCII output file
  CLOSE(ISCE)

  ! nullify pointers in the context setter
  call clear_sce_context()

  ! deallocate SCE and transformation arrays
  DEALLOCATE(APAR_MSP, BL_MSP, BU_MSP)
  DEALLOCATE(APAR_PHYS_MSP, BL_PHYS_MSP, BU_PHYS_MSP)
  DEALLOCATE(TRANSFORM_CODES)

  end subroutine sce_driver

end module sce_driver_MODULE
