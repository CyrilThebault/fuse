MODULE nr

  INTERFACE
    SUBROUTINE fdjac(x,fvec,df)
    USE nrtype
    REAL(WP), DIMENSION(:), INTENT(IN) :: fvec
    REAL(WP), DIMENSION(:), INTENT(INOUT) :: x
    REAL(WP), DIMENSION(:,:), INTENT(OUT) :: df
    END SUBROUTINE fdjac
  END INTERFACE

  INTERFACE gammln
    FUNCTION gammln_s(xx)
    USE nrtype
    REAL(WP), INTENT(IN) :: xx
    REAL(WP) :: gammln_s
    END FUNCTION gammln_s
    FUNCTION gammln_v(xx)
    USE nrtype
    REAL(WP), DIMENSION(:), INTENT(IN) :: xx
    REAL(WP), DIMENSION(size(xx)) :: gammln_v
    END FUNCTION gammln_v
  END INTERFACE

  INTERFACE gammp
    FUNCTION gammp_s(a,x)
    USE nrtype
    REAL(WP), INTENT(IN) :: a,x
    REAL(WP) :: gammp_s
    END FUNCTION gammp_s
    FUNCTION gammp_v(a,x)
    USE nrtype
    REAL(WP), DIMENSION(:), INTENT(IN) :: a,x
    REAL(WP), DIMENSION(size(a)) :: gammp_v
    END FUNCTION gammp_v
  END INTERFACE
  
  INTERFACE gammq
    FUNCTION gammq_s(a,x)
    USE nrtype
    REAL(WP), INTENT(IN) :: a,x
    REAL(WP) :: gammq_s
    END FUNCTION gammq_s
    FUNCTION gammq_v(a,x)
    USE nrtype
    REAL(WP), DIMENSION(:), INTENT(IN) :: a,x
    REAL(WP), DIMENSION(size(a)) :: gammq_v
    END FUNCTION gammq_v
  END INTERFACE
  
  INTERFACE gcf
    FUNCTION gcf_s(a,x,gln)
    USE nrtype
    REAL(WP), INTENT(IN) :: a,x
    REAL(WP), OPTIONAL, INTENT(OUT) :: gln
    REAL(WP) :: gcf_s
    END FUNCTION gcf_s
    FUNCTION gcf_v(a,x,gln)
    USE nrtype
    REAL(WP), DIMENSION(:), INTENT(IN) :: a,x
    REAL(WP), DIMENSION(:), OPTIONAL, INTENT(OUT) :: gln
    REAL(WP), DIMENSION(size(a)) :: gcf_v
    END FUNCTION gcf_v
  END INTERFACE
  
  INTERFACE gser
    FUNCTION gser_s(a,x,gln)
    USE nrtype
    REAL(WP), INTENT(IN) :: a,x
    REAL(WP), OPTIONAL, INTENT(OUT) :: gln
    REAL(WP) :: gser_s
    END FUNCTION gser_s
    FUNCTION gser_v(a,x,gln)
    USE nrtype
    REAL(WP), DIMENSION(:), INTENT(IN) :: a,x
    REAL(WP), DIMENSION(:), OPTIONAL, INTENT(OUT) :: gln
    REAL(WP), DIMENSION(size(a)) :: gser_v
    END FUNCTION gser_v
  END INTERFACE
  
  INTERFACE
    SUBROUTINE lnsrch(xold,fold,g,p,x,f,stpmax,check,func)
    USE nrtype
    REAL(WP), DIMENSION(:), INTENT(IN) :: xold,g
    REAL(WP), DIMENSION(:), INTENT(INOUT) :: p
    REAL(WP), INTENT(IN) :: fold,stpmax
    REAL(WP), DIMENSION(:), INTENT(OUT) :: x
    REAL(WP), INTENT(OUT) :: f
    LOGICAL(LGT), INTENT(OUT) :: check
    INTERFACE
      FUNCTION func(x)
      USE nrtype
      REAL(WP) :: func
      REAL(WP), DIMENSION(:), INTENT(IN) :: x
      END FUNCTION func
    END INTERFACE
    END SUBROUTINE lnsrch
  END INTERFACE

  INTERFACE
    SUBROUTINE lubksb(a,indx,b)
    USE nrtype
    REAL(WP), DIMENSION(:,:), INTENT(IN) :: a
    INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
    REAL(WP), DIMENSION(:), INTENT(INOUT) :: b
    END SUBROUTINE lubksb
  END INTERFACE
  
  INTERFACE
    SUBROUTINE ludcmp(a,indx,d)
    USE nrtype
    REAL(WP), DIMENSION(:,:), INTENT(INOUT) :: a
    INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: indx
    REAL(WP), INTENT(OUT) :: d
    END SUBROUTINE ludcmp
  END INTERFACE
  
  INTERFACE pythag
    FUNCTION pythag_sp(a,b)
    USE nrtype
    REAL(WP), INTENT(IN) :: a,b
    REAL(WP) :: pythag_sp
    END FUNCTION pythag_sp
  END INTERFACE
  
  INTERFACE svbksb
    SUBROUTINE svbksb_sp(u,w,v,b,x)
    USE nrtype
    REAL(WP), DIMENSION(:,:), INTENT(IN) :: u,v
    REAL(WP), DIMENSION(:), INTENT(IN) :: w,b
    REAL(WP), DIMENSION(:), INTENT(OUT) :: x
    END SUBROUTINE svbksb_sp
  END INTERFACE
  
  INTERFACE svdcmp
    SUBROUTINE svdcmp_sp(a,w,v)
    USE nrtype
    REAL(WP), DIMENSION(:,:), INTENT(INOUT) :: a
    REAL(WP), DIMENSION(:), INTENT(OUT) :: w
    REAL(WP), DIMENSION(:,:), INTENT(OUT) :: v
    END SUBROUTINE svdcmp_sp
  END INTERFACE

END MODULE nr
