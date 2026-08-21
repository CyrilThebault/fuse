  FUNCTION gcf_s(a,x,gln)
  USE nrtype; USE nrutil, ONLY : nrerror
  USE nr, ONLY : gammln
  IMPLICIT NONE
  REAL(WP), INTENT(IN) :: a,x
  REAL(WP), OPTIONAL, INTENT(OUT) :: gln
  REAL(WP) :: gcf_s
  INTEGER(I4B), PARAMETER :: ITMAX=100
  REAL(WP), PARAMETER :: EPS=epsilon(x),FPMIN=tiny(x)/EPS
  INTEGER(I4B) :: i
  REAL(WP) :: an,b,c,d,del,h
  if (x == 0.0_wp) then
    gcf_s=1._wp
    RETURN
  end if
  b=x+1.0_wp-a
  c=1.0_wp/FPMIN
  d=1.0_wp/b
  h=d
  do i=1,ITMAX
    an=-i*(i-a)
    b=b+2.0_wp
    d=an*d+b
    if (abs(d) < FPMIN) d=FPMIN
    c=b+an/c
    if (abs(c) < FPMIN) c=FPMIN
    d=1.0_wp/d
    del=d*c
    h=h*del
    if (abs(del-1.0_wp) <= EPS) exit
  end do
  if (i > ITMAX) call nrerror('a too large, ITMAX too small in gcf_s')
  if (present(gln)) then
    gln=gammln(a)
    gcf_s=exp(-x+a*log(x)-gln)*h
  else
    gcf_s=exp(-x+a*log(x)-gammln(a))*h
  end if
  END FUNCTION gcf_s


  FUNCTION gcf_v(a,x,gln)
  USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
  USE nr, ONLY : gammln
  IMPLICIT NONE
  REAL(WP), DIMENSION(:), INTENT(IN) :: a,x
  REAL(WP), DIMENSION(:), OPTIONAL, INTENT(OUT) :: gln
  REAL(WP), DIMENSION(size(a)) :: gcf_v
  INTEGER(I4B), PARAMETER :: ITMAX=100
  REAL(WP), PARAMETER :: EPS=epsilon(x),FPMIN=tiny(x)/EPS
  INTEGER(I4B) :: i
  REAL(WP), DIMENSION(size(a)) :: an,b,c,d,del,h
  LOGICAL(LGT), DIMENSION(size(a)) :: converged,zero
  i=assert_eq(size(a),size(x),'gcf_v')
  zero=(x == 0.0_wp)
  where (zero)
    gcf_v=1.0_wp
  elsewhere
    b=x+1.0_wp-a
    c=1.0_wp/FPMIN
    d=1.0_wp/b
    h=d
  end where
  converged=zero
  do i=1,ITMAX
    where (.not. converged)
      an=-i*(i-a)
      b=b+2.0_wp
      d=an*d+b
      d=merge(FPMIN,d, abs(d)<FPMIN )
      c=b+an/c
      c=merge(FPMIN,c, abs(c)<FPMIN )
      d=1.0_wp/d
      del=d*c
      h=h*del
      converged = (abs(del-1.0_wp)<=EPS)
    end where
    if (all(converged)) exit
  end do
  if (i > ITMAX) call nrerror('a too large, ITMAX too small in gcf_v')
  if (present(gln)) then
    if (size(gln) < size(a)) call &
      nrerror('gser: Not enough space for gln')
    gln=gammln(a)
    where (.not. zero) gcf_v=exp(-x+a*log(x)-gln)*h
  else
    where (.not. zero) gcf_v=exp(-x+a*log(x)-gammln(a))*h
  end if
  END FUNCTION gcf_v
