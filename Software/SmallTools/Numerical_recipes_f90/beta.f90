	FUNCTION beta_s(z,w)
	USE nrtype
	USE nr, ONLY : gammln
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: z,w
	REAL(SP) :: beta_s
	beta_s=exp(gammln(z)+gammln(w)-gammln(z+w))
	END FUNCTION beta_s


	FUNCTION beta_v(z,w)
	USE nrtype; USE nrutil, ONLY : assert_eq
	USE nr, ONLY : gammln
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: z,w
	REAL(SP), DIMENSION(size(z)) :: beta_v
	INTEGER(I4B) :: ndum
	ndum=assert_eq(size(z),size(w),'beta_v')
	beta_v=exp(gammln(z)+gammln(w)-gammln(z+w))
	END FUNCTION beta_v
