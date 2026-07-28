  FUNCTION gammp_s(a,x)
  USE nrtype; USE nrutil, ONLY : assert
  USE nr, ONLY : gcf,gser
  IMPLICIT NONE
  REAL(WP), INTENT(IN) :: a,x
  REAL(WP) :: gammp_s
  call assert( x >= 0.0,  a > 0.0, 'gammp_s args')
  if (x<a+1.0_wp) then
    gammp_s=gser(a,x)
  else
    gammp_s=1.0_wp-gcf(a,x)
  end if

  END FUNCTION gammp_s


  FUNCTION gammp_v(a,x)
  USE nrtype; USE nrutil, ONLY : assert,assert_eq
  USE nr, ONLY : gcf,gser
  IMPLICIT NONE
  REAL(WP), DIMENSION(:), INTENT(IN) :: a,x
  REAL(WP), DIMENSION(size(x)) :: gammp_v
  LOGICAL(LGT), DIMENSION(size(x)) :: mask
  INTEGER(I4B) :: ndum
  ndum=assert_eq(size(a),size(x),'gammp_v')
  call assert( all(x >= 0.0),  all(a > 0.0), 'gammp_v args')
  mask = (x<a+1.0_wp)
  gammp_v=merge(gser(a,merge(x,0.0_wp,mask)), &
    1.0_wp-gcf(a,merge(x,0.0_wp,.not. mask)),mask)
  END FUNCTION gammp_v
