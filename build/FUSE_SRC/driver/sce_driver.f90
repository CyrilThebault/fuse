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
  REAL(MSP), DIMENSION(:), ALLOCATABLE   :: URAND_MSP   ! vector of quasi-random numbers U[0,1]
  INTEGER(I4B)                           :: NOPT    ! number of parameters to be optimized
  INTEGER(I4B)                           :: NGS     ! # complexes in the initial population
  INTEGER(I4B)                           :: NPG     ! # points in each complex
  INTEGER(I4B)                           :: NPS     ! # points in a sub-complex
  INTEGER(I4B)                           :: NSPL    ! # evolution steps allowed for each complex before shuffling
  INTEGER(I4B)                           :: MINGS   ! minimum number of complexes required
  INTEGER(I4B)                           :: INIFLG  ! 1 = include initial point in the population
  INTEGER(I4B)                           :: IPRINT  ! 0 = supress printing
  INTEGER(I4B)                           :: ISCE    ! unit number for SCE write
  integer(i4b)                           :: NUMPSET ! number of parameter sets
  REAL(MSP)                              :: FUNCTN  ! function name for the model run
  INTEGER(KIND=4)                        :: ISEED   ! seed for the random sequence

  NOPT   =  NUMPAR         ! number of parameters to be optimized (NUMPAR in module multiparam)
  NGS    =     10          ! number of complexes in the initial population
  NPG    =  2*NOPT + 1     ! number of points in each complex
  NPS    =    NOPT + 1     ! number of points in a sub-complex
  NSPL   =  2*NOPT + 1     ! number of evolution steps allowed for each complex before shuffling
  MINGS  =  NGS            ! minimum number of complexes required
  INIFLG =  1              ! 1 = include initial point in the population
  IPRINT =  1              ! 0 = supress printing

  NUMPSET=1.2*MAXN         ! will be used to define the parameter set dimension of the NetCDF files
                           ! using 1.2MAXN since the final number of parameter sets produced by SCE is unknown

  ! convert from WP used in FUSE to MSP used in SCE
  ALLOCATE(APAR_MSP(NUMPAR), BL_MSP(NUMPAR), BU_MSP(NUMPAR))
  APAR_MSP=APAR; BL_MSP=BL; BU_MSP=BU

  ! pass the FUSE structures to the context setter
  ! NOTE: in sce_context_set, info/work/domain have the target attribute so can point to them
  call set_sce_context(info, work, domain)

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

  ! deallocate space for real32 vectors
  DEALLOCATE(APAR_MSP, BL_MSP, BU_MSP)
  
  end subroutine sce_driver

end module sce_driver_MODULE
