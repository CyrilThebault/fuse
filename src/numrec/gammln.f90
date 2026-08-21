  FUNCTION gammln_s(xx)
  USE nrtype; USE nrutil, ONLY : arth,assert
  IMPLICIT NONE
  REAL(WP), INTENT(IN) :: xx
  REAL(WP) :: gammln_s
  REAL(WP) :: tmp,x
  REAL(WP) :: stp = 2.5066282746310005_wp
  REAL(WP), DIMENSION(6) :: coef = (/76.18009172947146_wp,&
    -86.50532032941677_wp,24.01409824083091_wp,&
    -1.231739572450155_wp,0.1208650973866179e-2_wp,&
    -0.5395239384953e-5_wp/)
  call assert(xx > 0.0, 'gammln_s arg')
  x=xx
  tmp=x+5.5_wp
  tmp=(x+0.5_wp)*log(tmp)-tmp
  gammln_s=tmp+log(stp*(1.000000000190015_wp+&
    sum(coef(:)/arth(x+1.0_wp,1.0_wp,size(coef))))/x)
  END FUNCTION gammln_s


  FUNCTION gammln_v(xx)
  USE nrtype; USE nrutil, ONLY: assert
  IMPLICIT NONE
  INTEGER(I4B) :: i
  REAL(WP), DIMENSION(:), INTENT(IN) :: xx
  REAL(WP), DIMENSION(size(xx)) :: gammln_v
  REAL(WP), DIMENSION(size(xx)) :: ser,tmp,x,y
  REAL(WP) :: stp = 2.5066282746310005_wp
  REAL(WP), DIMENSION(6) :: coef = (/76.18009172947146_wp,&
    -86.50532032941677_wp,24.01409824083091_wp,&
    -1.231739572450155_wp,0.1208650973866179e-2_wp,&
    -0.5395239384953e-5_wp/)
  if (size(xx) == 0) RETURN
  call assert(all(xx > 0.0), 'gammln_v arg')
  x=xx
  tmp=x+5.5_wp
  tmp=(x+0.5_wp)*log(tmp)-tmp
  ser=1.000000000190015_wp
  y=x
  do i=1,size(coef)
    y=y+1.0_wp
    ser=ser+coef(i)/y
  end do
  gammln_v=tmp+log(stp*ser/x)
  END FUNCTION gammln_v
