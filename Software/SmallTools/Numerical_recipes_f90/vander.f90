	FUNCTION vander(x,q)
	USE nrtype; USE nrutil, ONLY : assert_eq,outerdiff
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), INTENT(IN) :: x,q
	REAL(DP), DIMENSION(size(x)) :: vander
	REAL(DP), DIMENSION(size(x)) :: c
	REAL(DP), DIMENSION(size(x),size(x)) :: a
	INTEGER(I4B) :: i,n
	n=assert_eq(size(x),size(q),'vander')
	if (n == 1) then
		vander(1)=q(1)
	else
		c(:)=0.0
		c(n)=-x(1)
		do i=2,n
			c(n+1-i:n-1)=c(n+1-i:n-1)-x(i)*c(n+2-i:n)
			c(n)=c(n)-x(i)
		end do
		a(:,:)=outerdiff(x,x)
		vander(:)=product(a,dim=2,mask=(a /= 0.0))
		a(:,1)=-c(1)/x(:)
		do i=2,n
			a(:,i)=-(c(i)-a(:,i-1))/x(:)
		end do
		vander(:)=matmul(a,q)/vander(:)
	end if
	END FUNCTION vander
