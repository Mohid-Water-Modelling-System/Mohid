	SUBROUTINE rsolv(a,d,b)
	USE nrtype; USE nrutil, ONLY : assert_eq
	IMPLICIT NONE
	REAL(SP), DIMENSION(:,:), INTENT(IN) :: a
	REAL(SP), DIMENSION(:), INTENT(IN) :: d
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: b
	INTEGER(I4B) :: i,n
	n=assert_eq(size(a,1),size(a,2),size(b),size(d),'rsolv')
	b(n)=b(n)/d(n)
	do i=n-1,1,-1
		b(i)=(b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/d(i)
	end do
	END SUBROUTINE rsolv
