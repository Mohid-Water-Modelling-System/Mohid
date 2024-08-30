	FUNCTION fpoly(x,n)
	USE nrtype; USE nrutil, ONLY : geop
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	INTEGER(I4B), INTENT(IN) :: n
	REAL(SP), DIMENSION(n) :: fpoly
	fpoly=geop(1.0_sp,x,n)
	END FUNCTION fpoly
