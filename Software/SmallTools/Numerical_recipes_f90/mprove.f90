	SUBROUTINE mprove(a,alud,indx,b,x)
	USE nrtype; USE nrutil, ONLY : assert_eq
	USE nr, ONLY : lubksb
	IMPLICIT NONE
	REAL(SP), DIMENSION(:,:), INTENT(IN) :: a,alud
	INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
	REAL(SP), DIMENSION(:), INTENT(IN) :: b
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: x
	INTEGER(I4B) :: ndum
	REAL(SP), DIMENSION(size(a,1)) :: r
	ndum=assert_eq((/size(a,1),size(a,2),size(alud,1),size(alud,2),size(b),&
		size(x),size(indx)/),'mprove')
	r=matmul(real(a,dp),real(x,dp))-real(b,dp)
	call lubksb(alud,indx,r)
	x=x-r
	END SUBROUTINE mprove
