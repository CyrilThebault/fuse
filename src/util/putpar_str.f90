module putpar_str_module

  implicit none
  private

  public :: putpar_str

contains

  subroutine putpar_str(metadat, parname, parmeta, ierr, message)

    ! variable types
    use nrtype, only: wp, i4b
    use work_types, only: fuse_work

    ! legacy parameter metadata structures
    use multiparam_types, only: paratt, parinfo

    implicit none

    ! input
    type(paratt),  intent(in)    :: metadat
    character(*),  intent(in)    :: parname
    type(parinfo), intent(inout) :: parmeta

    ! error control
    integer(i4b)           , intent(out)   :: ierr
    character(*)           , intent(out)   :: message

    ierr    = 0
    message = "putpar_str/"

    SELECTCASE(TRIM(PARNAME))
      CASE('RFERR_ADD');  PARMETA%RFERR_ADD = METADAT 
      CASE('RFERR_MLT');  PARMETA%RFERR_MLT = METADAT
      CASE('RFH1_MEAN');  PARMETA%RFH1_MEAN = METADAT
      CASE('RFH2_SDEV');  PARMETA%RFH2_SDEV = METADAT
      CASE('RH1P_MEAN');  PARMETA%RH1P_MEAN = METADAT
      CASE('RH1P_SDEV');  PARMETA%RH1P_SDEV = METADAT
      CASE('RH2P_MEAN');  PARMETA%RH2P_MEAN = METADAT
      CASE('RH2P_SDEV');  PARMETA%RH2P_SDEV = METADAT
      CASE('REFSINT_0');  PARMETA%REFSINT_0 = METADAT
      CASE('MAXWATR_1');  PARMETA%MAXWATR_1 = METADAT
      CASE('MAXWATR_2');  PARMETA%MAXWATR_2 = METADAT
      CASE('FRACTEN');    PARMETA%FRACTEN   = METADAT
      CASE('FRCHZNE');    PARMETA%FRCHZNE   = METADAT
      CASE('FPRIMQB');    PARMETA%FPRIMQB   = METADAT
      CASE('RTFRAC1');    PARMETA%RTFRAC1   = METADAT
      CASE('PERCRTE');    PARMETA%PERCRTE   = METADAT
      CASE('PERCEXP');    PARMETA%PERCEXP   = METADAT
      CASE('SACPMLT');    PARMETA%SACPMLT   = METADAT
      CASE('SACPEXP');    PARMETA%SACPEXP   = METADAT
      CASE('PERCFRAC');   PARMETA%PERCFRAC  = METADAT
      CASE('FRACLOWZ');   PARMETA%FRACLOWZ  = METADAT
      CASE('IFLWRTE');    PARMETA%IFLWRTE   = METADAT
      CASE('BASERTE');    PARMETA%BASERTE   = METADAT
      CASE('QB_POWR');    PARMETA%QB_POWR   = METADAT
      CASE('QB_PRMS');    PARMETA%QB_PRMS   = METADAT
      CASE('QBRATE_2A');  PARMETA%QBRATE_2A = METADAT
      CASE('QBRATE_2B');  PARMETA%QBRATE_2B = METADAT
      CASE('SAREAMAX');   PARMETA%SAREAMAX  = METADAT
      CASE('AXV_BEXP');   PARMETA%AXV_BEXP  = METADAT
      CASE('LOGLAMB');    PARMETA%LOGLAMB   = METADAT
      CASE('TISHAPE');    PARMETA%TISHAPE   = METADAT
      CASE('TIMEDELAY');  PARMETA%TIMEDELAY = METADAT
      CASE('MBASE');      PARMETA%MBASE     = METADAT
      CASE('MFMAX');      PARMETA%MFMAX     = METADAT
      CASE('MFMIN');      PARMETA%MFMIN     = METADAT
      CASE('PXTEMP');     PARMETA%PXTEMP    = METADAT
      CASE('OPG');        PARMETA%OPG       = METADAT
      CASE('LAPSE');      PARMETA%LAPSE     = METADAT
      CASE DEFAULT
        message = trim(message)//'parameter name "'//trim(parname)//'" does not exist'
        ierr=10; return
    ENDSELECT

  end subroutine putpar_str

end module putpar_str_module
