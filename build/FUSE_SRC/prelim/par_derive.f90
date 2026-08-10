module PAR_DERIVE_module

  use work_types, only: fuse_param

  implicit none
  
  private
  public :: PAR_DERIVE

contains

  SUBROUTINE PAR_DERIVE(parStruct)
  ! ---------------------------------------------------------------------------------------
  ! Creator:
  ! --------
  ! Martyn Clark, 2007
  ! ---------------------------------------------------------------------------------------
  ! Purpose:
  ! --------
  ! Computes derived model parameters (bucket sizes, etc.)
  ! ---------------------------------------------------------------------------------------
  ! Modules Modified:
  ! -----------------
  ! MODULE multiparam -- model parameters stored in MODULE multiparam
  ! ---------------------------------------------------------------------------------------
  
  ! model definition structures
  USE model_defn, ONLY: SMODL
  USE model_defnames

  ! shared data
  USE multiparam, ONLY: MPARAM, DPARAM  ! model parameter structures
  
  IMPLICIT NONE

  type(fuse_param), intent(inout)  :: parStruct
  
  CALL BUCKETSIZE()        ! compute bucket size
  CALL MEAN_TIPOW()        ! mean of the power-transformed topo index
  CALL QBSATURATN()        ! compute baseflow at saturation (used in the SAC percolation model)
  
  IF (SMODL%iESOIL.EQ.iopt_rootweight) DPARAM%RTFRAC2 = 1._WP - MPARAM%RTFRAC1
  
  parStruct%param_derive = DPARAM
  
  END SUBROUTINE PAR_DERIVE

end module PAR_DERIVE_module
