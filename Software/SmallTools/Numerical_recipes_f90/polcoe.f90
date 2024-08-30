	FUNCTION polcoe(x,y)
	USE nrtype; USE nrutil, ONLY : assert_eq,outerdiff
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
	REAL(SP), DIMENSION(size(x)) :: polcoe
	INTEGER(I4B) :: i,k,n
	REAL(SP), DIMENSION(size(x)) :: s
	REAL(SP), DIMENSION(size(x),size(x)) :: a
	n=assert_eq(size(x),size(y),'polcoe')
	s=0.0
	s(n)=-x(1)
	do i=2,n
		s(n+1-i:n-1)=s(n+1-i:n-1)-x(i)*s(n+2-i:n)
		s(n)=s(n)-x(i)
	end do
	a=outerdiff(x,x)
	polcoe=product(a,dim=2,mask=a /= 0.0)
	a(:,1)=-s(1)/x(:)
	do k=2,n
		a(:,k)=-(s(k)-a(:,k-1))/x(:)
	end do
	s=y/polcoe
	polcoe=matmul(s,a)
	END FUNCTION polcoe
