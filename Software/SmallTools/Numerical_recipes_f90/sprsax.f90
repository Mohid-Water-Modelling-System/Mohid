	SUBROUTINE sprsax_sp(sa,x,b)
	USE nrtype; USE nrutil, ONLY : assert_eq,scatter_add
	IMPLICIT NONE
	TYPE(sprs2_sp), INTENT(IN) :: sa
	REAL(SP), DIMENSION (:), INTENT(IN) :: x
	REAL(SP), DIMENSION (:), INTENT(OUT) :: b
	INTEGER(I4B) :: ndum
	ndum=assert_eq(sa%n,size(x),size(b),'sprsax_sp')
	b=0.0_sp
	call scatter_add(b,sa%val*x(sa%jcol),sa%irow)
	END SUBROUTINE sprsax_sp

	SUBROUTINE sprsax_dp(sa,x,b)
	USE nrtype; USE nrutil, ONLY : assert_eq,scatter_add
	IMPLICIT NONE
	TYPE(sprs2_dp), INTENT(IN) :: sa
	REAL(DP), DIMENSION (:), INTENT(IN) :: x
	REAL(DP), DIMENSION (:), INTENT(OUT) :: b
	INTEGER(I4B) :: ndum
	ndum=assert_eq(sa%n,size(x),size(b),'sprsax_dp')
	b=0.0_dp
	call scatter_add(b,sa%val*x(sa%jcol),sa%irow)
	END SUBROUTINE sprsax_dp
