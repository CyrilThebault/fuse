module varextract_module

  use nrtype
  use iso_fortran_env, only: real32
  use domain_types, only: domain_data
  use fuse_globaldata, only: NA_VALUE_SP

  implicit none
  private
  public :: varextract_3d

contains

  subroutine varextract_3d(domain, varname, nspat1, nspat2, numtim, xout)
  ! ---------------------------------------------------------------------------------------
  USE model_numerix
  USE multiforce, only: aValid                                ! model validation data
  USE multistate, only: gState_3d                             ! model states
  USE multi_flux, only: w_flux_3d                             ! model fluxes
  USE multiroute, only: aroute_3d                             ! routed runoff
  implicit none

  type(domain_data), intent(in) :: domain
  character(*),      intent(in) :: varname
  integer(i4b),      intent(in) :: nspat1, nspat2, numtim

  ! NetCDF output buffer (matches NF90_FLOAT)
  real(real32), intent(out) :: xout(nspat1, nspat2, numtim)

  ! ---------------------------------------------------------------------------------------
  ! the length of the temporal dimension of the state variables (gState_3d and MBANDS_VAR_4d)
  ! is greater by one time step, so only keeping first numtim time steps, i.e. not writing
  ! last value the output file
  
  SELECT CASE (TRIM(VARNAME))
  
   ! extract forcing data
   CASE ('ppt')        ; xout = real(domain%force(:,:,1:numtim)%PPT ,        kind=real32)
   CASE ('temp')       ; xout = real(domain%force(:,:,1:numtim)%TEMP,        kind=real32) 
   CASE ('pet')        ; xout = real(domain%force(:,:,1:numtim)%PET ,        kind=real32) 
   
   ! extract response data
   ! TODO: Check this -- it is weird that obs q is 3d
   CASE ('obsq')       ; xout = real(aValid(:,:,1:numtim)%OBSQ,           kind=real32)
   
   ! extract model states
   CASE ('tens_1')     ; xout = real(gState_3d(:,:,1:numtim)%TENS_1 ,     kind=real32)
   CASE ('tens_1a')    ; xout = real(gState_3d(:,:,1:numtim)%TENS_1A,     kind=real32) 
   CASE ('tens_1b')    ; xout = real(gState_3d(:,:,1:numtim)%TENS_1B,     kind=real32) 
   CASE ('free_1')     ; xout = real(gState_3d(:,:,1:numtim)%FREE_1 ,     kind=real32)
   CASE ('watr_1')     ; xout = real(gState_3d(:,:,1:numtim)%WATR_1 ,     kind=real32)
   CASE ('tens_2')     ; xout = real(gState_3d(:,:,1:numtim)%TENS_2 ,     kind=real32)
   CASE ('free_2')     ; xout = real(gState_3d(:,:,1:numtim)%FREE_2 ,     kind=real32)
   CASE ('free_2a')    ; xout = real(gState_3d(:,:,1:numtim)%FREE_2A,     kind=real32) 
   CASE ('free_2b')    ; xout = real(gState_3d(:,:,1:numtim)%FREE_2B,     kind=real32) 
   CASE ('watr_2')     ; xout = real(gState_3d(:,:,1:numtim)%WATR_2 ,     kind=real32)
   CASE ('swe_tot')    ; xout = real(gState_3d(:,:,1:numtim)%swe_tot,     kind=real32) 
   
   ! extract model fluxes
   CASE ('eff_ppt')    ; xout = real(W_FLUX_3d(:,:,1:numtim)%EFF_PPT    , kind=real32)
   CASE ('satarea')    ; xout = real(W_FLUX_3d(:,:,1:numtim)%SATAREA    , kind=real32)
   CASE ('qsurf')      ; xout = real(W_FLUX_3d(:,:,1:numtim)%QSURF      , kind=real32)
   CASE ('evap_1a')    ; xout = real(W_FLUX_3d(:,:,1:numtim)%EVAP_1A    , kind=real32)
   CASE ('evap_1b')    ; xout = real(W_FLUX_3d(:,:,1:numtim)%EVAP_1B    , kind=real32)
   CASE ('evap_1')     ; xout = real(W_FLUX_3d(:,:,1:numtim)%EVAP_1     , kind=real32)
   CASE ('evap_2')     ; xout = real(W_FLUX_3d(:,:,1:numtim)%EVAP_2     , kind=real32)
   CASE ('rchr2excs')  ; xout = real(W_FLUX_3d(:,:,1:numtim)%RCHR2EXCS  , kind=real32)
   CASE ('tens2free_1'); xout = real(W_FLUX_3d(:,:,1:numtim)%TENS2FREE_1, kind=real32) 
   CASE ('oflow_1')    ; xout = real(W_FLUX_3d(:,:,1:numtim)%OFLOW_1    , kind=real32)
   CASE ('tens2free_2'); xout = real(W_FLUX_3d(:,:,1:numtim)%TENS2FREE_2, kind=real32)  
   CASE ('qintf_1')    ; xout = real(W_FLUX_3d(:,:,1:numtim)%QINTF_1    , kind=real32)
   CASE ('qperc_12')   ; xout = real(W_FLUX_3d(:,:,1:numtim)%QPERC_12   , kind=real32)
   CASE ('qbase_2')    ; xout = real(W_FLUX_3d(:,:,1:numtim)%QBASE_2    , kind=real32)
   CASE ('qbase_2a')   ; xout = real(W_FLUX_3d(:,:,1:numtim)%QBASE_2A   , kind=real32)
   CASE ('qbase_2b')   ; xout = real(W_FLUX_3d(:,:,1:numtim)%QBASE_2B   , kind=real32)
   CASE ('oflow_2')    ; xout = real(W_FLUX_3d(:,:,1:numtim)%OFLOW_2    , kind=real32)
   CASE ('oflow_2a')   ; xout = real(W_FLUX_3d(:,:,1:numtim)%OFLOW_2A   , kind=real32)
   CASE ('oflow_2b')   ; xout = real(W_FLUX_3d(:,:,1:numtim)%OFLOW_2B   , kind=real32)
   
   ! extract extrapolation errors
   CASE ('err_tens_1') ; xout = real(W_FLUX_3d(:,:,1:numtim)%ERR_TENS_1 , kind=real32)
   CASE ('err_tens_1a'); xout = real(W_FLUX_3d(:,:,1:numtim)%ERR_TENS_1A, kind=real32)
   CASE ('err_tens_1b'); xout = real(W_FLUX_3d(:,:,1:numtim)%ERR_TENS_1B, kind=real32)
   CASE ('err_free_1') ; xout = real(W_FLUX_3d(:,:,1:numtim)%ERR_FREE_1 , kind=real32)
   CASE ('err_watr_1') ; xout = real(W_FLUX_3d(:,:,1:numtim)%ERR_WATR_1 , kind=real32)
   CASE ('err_tens_2') ; xout = real(W_FLUX_3d(:,:,1:numtim)%ERR_TENS_2 , kind=real32)
   CASE ('err_free_2') ; xout = real(W_FLUX_3d(:,:,1:numtim)%ERR_FREE_2 , kind=real32)
   CASE ('err_free_2a'); xout = real(W_FLUX_3d(:,:,1:numtim)%ERR_FREE_2A, kind=real32)
   CASE ('err_free_2b'); xout = real(W_FLUX_3d(:,:,1:numtim)%ERR_FREE_2B, kind=real32)
   CASE ('err_watr_2') ; xout = real(W_FLUX_3d(:,:,1:numtim)%ERR_WATR_2 , kind=real32)
   
   ! time check
   CASE ('chk_time')   ; xout = real(W_FLUX_3d(:,:,1:numtim)%CHK_TIME,    kind=real32)
   
   ! extract model runoff
   CASE ('q_instnt')   ; xout = real(AROUTE_3d(:,:,1:numtim)%Q_INSTNT,    kind=real32)
   CASE ('q_routed')   ; xout = real(AROUTE_3d(:,:,1:numtim)%Q_ROUTED,    kind=real32)
   
   ! extract information on numerical solution (shared in MODULE model_numerix)
   ! TODO: Check the need for this -- broadcasting scalars to the 3-d field
   CASE ('num_funcs')  ; xout = real(NUM_FUNCS    , kind=real32)
   CASE ('numjacobian'); xout = real(NUM_JACOBIAN , kind=real32)
   CASE ('sub_accept') ; xout = real(NUMSUB_ACCEPT, kind=real32)
   CASE ('sub_reject') ; xout = real(NUMSUB_REJECT, kind=real32)
   CASE ('sub_noconv') ; xout = real(NUMSUB_NOCONV, kind=real32)
   CASE ('max_iterns') ; xout = real(MAXNUM_ITERNS, kind=real32)

   ! default
   case default;         xout = NA_VALUE_SP 

  END SELECT
  
  ! ---------------------------------------------------------------------------------------
  END SUBROUTINE VAREXTRACT_3d

END MODULE VAREXTRACT_MODULE
