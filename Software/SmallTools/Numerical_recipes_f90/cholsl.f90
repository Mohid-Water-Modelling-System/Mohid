	SUBROUTINE cholsl(a,p,b,x)
	USE nrtype; USE nrutil, ONLY : assert_eq
	IMPLICIT NONE
	REAL(SP), DIMENSION(:,:), INTENT(IN) :: a
	REAL(SP), DIMENSION(:), INTENT(IN) :: p,b
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: x
	INTEGER(I4B) :: i,n
	n=assert_eq((/size(a,1),size(a,2),size(p),size(b),size(x)/),'cholsl')
	do i=1,n
		x(i)=(b(i)-dot_product(a(i,1:i-1),x(1:i-1)))/p(i)
	end do
	do i=n,1,-1
		x(i)=(x(i)-dot_product(a(i+1:n,i),x(i+1:n)))/p(i)
	end do
	END SUBROUTINE cholsl
