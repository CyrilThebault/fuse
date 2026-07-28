  FUNCTION pythag_sp(a,b)
  USE nrtype
  IMPLICIT NONE
  REAL(WP), INTENT(IN) :: a,b
  REAL(WP) :: pythag_sp
  REAL(WP) :: absa,absb
  absa=abs(a)
  absb=abs(b)
  if (absa > absb) then
    pythag_sp=absa*sqrt(1.0_wp+(absb/absa)**2)
  else
    if (absb == 0.0) then
      pythag_sp=0.0
    else
      pythag_sp=absb*sqrt(1.0_wp+(absa/absb)**2)
    end if
  end if
  END FUNCTION pythag_sp
