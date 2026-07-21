module get_bundle_module
  use nrtype
  use work_types, only: fuse_work
  USE model_defn, ONLY: NSTATE   ! TODO: update to new structures
  USE multiparam, ONLY: NUMPAR   ! TODO: update to new structures
  implicit none

contains

  subroutine get_bundle(fuseStruct)
  use multiforce, only: timDat
  use multiforce, only: mForce
  use multistate, only: mState
  use multi_flux, only: m_flux
  use multiparam, only: parMeta,mParam,dParam
  implicit none
  type(fuse_work), intent(inout) :: fuseStruct
  integer(i4b)                   :: iState
  integer(i4b)                   :: iParam

  ! populate fuse work structures
  fuseStruct%step%time         = timdat
  fuseStruct%step%force        = mForce
  fuseStruct%step%state0       = mState
  fuseStruct%step%state1       = mState
  fuseStruct%step%flux         = m_flux  ! initialized at zero

  fuseStruct%par%param_meta   = parMeta
  fuseStruct%par%param_adjust = mParam
  fuseStruct%par%param_derive = dParam

  ! initialize flux derivatives
  do iState=1,nState
   fuseStruct%adj%df_dS(iState) = m_flux ! initialized at zero
  end do

  ! initialize parameter derivatives
  do iParam=1,NUMPAR
   fuseStruct%adj%df_dPar(iParam) = m_flux ! initialized at zero
  end do

  end subroutine get_bundle


end module get_bundle_module
