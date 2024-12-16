	SUBROUTINE choldc(a,p)
	USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
	IMPLICIT NONE
	REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
	REAL(SP), DIMENSION(:), INTENT(OUT) :: p
	INTEGER(I4B) :: i,n
	REAL(SP) :: summ
	n=assert_eq(size(a,1),size(a,2),size(p),'choldc')
	do i=1,n
		summ=a(i,i)-dot_product(a(i,1:i-1),a(i,1:i-1))
		if (summ <= 0.0) call nrerror('choldc failed')
		p(i)=sqrt(summ)
		a(i+1:n,i)=(a(i,i+1:n)-matmul(a(i+1:n,1:i-1),a(i,1:i-1)))/p(i)
	end do
	END SUBROUTINE choldc
