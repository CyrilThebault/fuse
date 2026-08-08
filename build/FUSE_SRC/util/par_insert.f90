! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark
! Modified by Brian Henn to include snow model, 6/2013
! ---------------------------------------------------------------------------------------
module par_insert_module

  use nrtype,           only: i4b, wp
  use work_types,       only: fuse_param
  use multiparam_types, only: par_id

  implicit none
  private

  public :: put_parset
  public :: par_insert

contains

  !---------------------------------------------------------------------
  ! Insert a complete parameter vector into the FUSE parameter
  ! structures using the corresponding parameter names.
  !---------------------------------------------------------------------
  subroutine put_parset(parset, parnames, parStruct, ierr, message)

    implicit none

    real(wp),         intent(in)    :: parset(:)
    type(par_id),     intent(in)    :: parnames(:)
    type(fuse_param), intent(inout) :: parStruct
    
    integer(i4b),      intent(out)  :: ierr
    character(*),      intent(out)  :: message

    integer(i4b)                    :: ipar
    character(len=256)              :: cmessage         ! error message of downwind routine


    ierr=0
    message='put_parset/'

    if (size(parset) /= size(parnames)) then
      message=trim(message)//'inconsistent parameter dimensions'
      ierr=20; return
    end if

    do ipar = 1, size(parset)

      call par_insert(               &
          parset(ipar),              & 
          parnames(ipar)%parname,    &
          parStruct,                 &
          ierr, cmessage)

      if(ierr/=0)then
        message=trim(message)//trim(cmessage)
        return
      endif

    end do
  end subroutine put_parset

  !---------------------------------------------------------------------
  ! Insert one parameter value into the FUSE parameter structures.
  !---------------------------------------------------------------------
  subroutine par_insert(xvar, parname, parStruct, ierr, message)

    implicit none

    real(wp),         intent(in)    :: xvar
    character(*),     intent(in)    :: parname
    type(fuse_param), intent(inout) :: parStruct

    integer(i4b),      intent(out)  :: ierr
    character(*),      intent(out)  :: message

    ierr=0
    message='put_insert/'

    SELECTCASE(TRIM(PARNAME))

    ! adjustable model parameters
    CASE('RFERR_ADD');  parStruct%param_adjust%RFERR_ADD  = XVAR
    CASE('RFERR_MLT');  parStruct%param_adjust%RFERR_MLT  = XVAR
    CASE('RFH1_MEAN');  parStruct%param_adjust%RFH1_MEAN  = XVAR
    CASE('RFH2_SDEV');  parStruct%param_adjust%RFH2_SDEV  = XVAR
    CASE('RH1P_MEAN');  parStruct%param_adjust%RH1P_MEAN  = XVAR
    CASE('RH1P_SDEV');  parStruct%param_adjust%RH1P_SDEV  = XVAR
    CASE('RH2P_MEAN');  parStruct%param_adjust%RH2P_MEAN  = XVAR
    CASE('RH2P_SDEV');  parStruct%param_adjust%RH2P_SDEV  = XVAR
    CASE('REFSINT_0');  parStruct%param_adjust%REFSINT_0  = XVAR
    CASE('MAXWATR_1');  parStruct%param_adjust%MAXWATR_1  = XVAR
    CASE('MAXWATR_2');  parStruct%param_adjust%MAXWATR_2  = XVAR
    CASE('FRACTEN');    parStruct%param_adjust%FRACTEN    = XVAR
    CASE('FRCHZNE');    parStruct%param_adjust%FRCHZNE    = XVAR
    CASE('FPRIMQB');    parStruct%param_adjust%FPRIMQB    = XVAR
    CASE('RTFRAC1');    parStruct%param_adjust%RTFRAC1    = XVAR
    CASE('PERCRTE');    parStruct%param_adjust%PERCRTE    = XVAR
    CASE('PERCEXP');    parStruct%param_adjust%PERCEXP    = XVAR
    CASE('SACPMLT');    parStruct%param_adjust%SACPMLT    = XVAR
    CASE('SACPEXP');    parStruct%param_adjust%SACPEXP    = XVAR
    CASE('PERCFRAC');   parStruct%param_adjust%PERCFRAC   = XVAR
    CASE('FRACLOWZ');   parStruct%param_adjust%FRACLOWZ   = XVAR
    CASE('IFLWRTE');    parStruct%param_adjust%IFLWRTE    = XVAR
    CASE('BASERTE');    parStruct%param_adjust%BASERTE    = XVAR
    CASE('QB_POWR');    parStruct%param_adjust%QB_POWR    = XVAR
    CASE('QB_PRMS');    parStruct%param_adjust%QB_PRMS    = XVAR
    CASE('QBRATE_2A');  parStruct%param_adjust%QBRATE_2A  = XVAR
    CASE('QBRATE_2B');  parStruct%param_adjust%QBRATE_2B  = XVAR
    CASE('SAREAMAX');   parStruct%param_adjust%SAREAMAX   = XVAR
    CASE('AXV_BEXP');   parStruct%param_adjust%AXV_BEXP   = XVAR
    CASE('LOGLAMB');    parStruct%param_adjust%LOGLAMB    = XVAR
    CASE('TISHAPE');    parStruct%param_adjust%TISHAPE    = XVAR
    CASE('TIMEDELAY');  parStruct%param_adjust%TIMEDELAY  = XVAR
    CASE('MBASE');      parStruct%param_adjust%MBASE      = XVAR
    CASE('MFMAX');      parStruct%param_adjust%MFMAX      = XVAR
    CASE('MFMIN');      parStruct%param_adjust%MFMIN      = XVAR
    CASE('PXTEMP');     parStruct%param_adjust%PXTEMP     = XVAR
    CASE('OPG');        parStruct%param_adjust%OPG        = XVAR
    CASE('LAPSE');      parStruct%param_adjust%LAPSE      = XVAR

      ! derived parameters
    CASE('MAXTENS_1');  parStruct%param_derive%MAXTENS_1  = XVAR
    CASE('MAXTENS_1A'); parStruct%param_derive%MAXTENS_1A = XVAR
    CASE('MAXTENS_1B'); parStruct%param_derive%MAXTENS_1B = XVAR
    CASE('MAXFREE_1');  parStruct%param_derive%MAXFREE_1  = XVAR
    CASE('MAXTENS_2');  parStruct%param_derive%MAXTENS_2  = XVAR
    CASE('MAXFREE_2');  parStruct%param_derive%MAXFREE_2  = XVAR
    CASE('MAXFREE_2A'); parStruct%param_derive%MAXFREE_2A = XVAR
    CASE('MAXFREE_2B'); parStruct%param_derive%MAXFREE_2B = XVAR
    CASE('QBSAT');      parStruct%param_derive%QBSAT      = XVAR
    CASE('RTFRAC2');    parStruct%param_derive%RTFRAC2    = XVAR
    CASE('POWLAMB');    parStruct%param_derive%POWLAMB    = XVAR
    CASE('MAXPOW');     parStruct%param_derive%MAXPOW     = XVAR

    CASE DEFAULT
      message=trim(message)//'parameter name does not exist '
      ierr=20; return
    
    ENDSELECT
  
  end subroutine par_insert

end module par_insert_module
