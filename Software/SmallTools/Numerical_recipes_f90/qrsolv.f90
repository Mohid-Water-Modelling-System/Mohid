	SUBROUTINE qrsolv(a,c,d,b)
	USE nrtype; USE nrutil, ONLY : assert_eq
	USE nr, ONLY : rsolv
	IMPLICIT NONE
	REAL(SP), DIMENSION(:,:), INTENT(IN) :: a
	REAL(SP), DIMENSION(:), INTENT(IN) :: c,d
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: b
	INTEGER(I4B) :: j,n
	REAL(SP) :: tau
	n=assert_eq((/size(a,1),size(a,2),size(b),size(c),size(d)/),'qrsolv')
	do j=1,n-1
		tau=dot_product(a(j:n,j),b(j:n))/c(j)
		b(j:n)=b(j:n)-tau*a(j:n,j)
	end do
	call rsolv(a,d,b)
	END SUBROUTINE qrsolv
