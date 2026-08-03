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

  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: TRANSFORM_CODES

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

  call setup_parameter_transforms( APAR,     BL,     BU,       &
                                   APAR_MSP, BL_MSP, BU_MSP,   &
                                   TRANSFORM_CODES,            &
                                   IERR, MESSAGE )

  if (IERR /= 0) then
    write(*,'(A)') trim(MESSAGE)
    stop 'Unable to set up parameter transformations'
  end if

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
  DEALLOCATE(TRANSFORM_CODES)

  end subroutine sce_driver

  ! ------------------------------------------------------------------------
  ! Prepare parameter transformations for SCE optimization.
  ! ------------------------------------------------------------------------
  subroutine setup_parameter_transforms( apar_phys,   bl_phys,   bu_phys,    &
                                         apar_search, bl_search, bu_search,  &
                                         transform_codes, ierr, message)

    use multiparam, only: NUMPAR, LPARAM, PARATT
    use getpar_str_module, only: getpar_str
    use parameter_transform_module, only: validate_transform
    use parameter_transform_module, only: vector_to_search_space

    implicit none

    real(wp), intent(in) :: apar_phys(:)
    real(wp), intent(in) :: bl_phys(:)
    real(wp), intent(in) :: bu_phys(:)

    real(MSP), allocatable, intent(out) :: apar_search(:)
    real(MSP), allocatable, intent(out) :: bl_search(:)
    real(MSP), allocatable, intent(out) :: bu_search(:)

    integer(I4B), allocatable, intent(out) :: transform_codes(:)

    integer(I4B), intent(out) :: ierr
    character(len=*), intent(out) :: message

    real(MSP), allocatable :: apar_phys_msp(:)
    real(MSP), allocatable :: bl_phys_msp(:)
    real(MSP), allocatable :: bu_phys_msp(:)

    type(PARATT) :: param_meta
    integer(I4B) :: ipar

    ierr = 0
    message = ''

    allocate(apar_phys_msp(NUMPAR))
    allocate(bl_phys_msp(NUMPAR))
    allocate(bu_phys_msp(NUMPAR))

    allocate(apar_search(NUMPAR))
    allocate(bl_search(NUMPAR))
    allocate(bu_search(NUMPAR))

    allocate(transform_codes(NUMPAR))

    apar_phys_msp = apar_phys
    bl_phys_msp   = bl_phys
    bu_phys_msp   = bu_phys

    do ipar = 1, NUMPAR

      call getpar_str(LPARAM(ipar)%PARNAME, param_meta)

      transform_codes(ipar) = param_meta%PARVTN

      call validate_transform(LPARAM(ipar)%PARNAME, bl_phys_msp(ipar), &
                              apar_phys_msp(ipar), bu_phys_msp(ipar),  &
                              transform_codes(ipar), ierr, message)

      if (ierr /= 0) then
        message = "Parameter validation: "//trim(message)
        return
      end if

    end do

    call vector_to_search_space(apar_phys_msp, transform_codes, apar_search, ierr, message)

    if (ierr /= 0) then
      message = 'Initial parameter set: '//trim(message)
      return
    end if

    call vector_to_search_space(bl_phys_msp, transform_codes, bl_search, ierr, message)

    if (ierr /= 0) then
      message = 'Lower parameter bounds: '//trim(message)
      return
    end if

    call vector_to_search_space(bu_phys_msp, transform_codes, bu_search, ierr, message)

    if (ierr /= 0) then
      message = 'Upper parameter bounds: '//trim(message)
      return
    end if

  end subroutine setup_parameter_transforms

end module sce_driver_MODULE
