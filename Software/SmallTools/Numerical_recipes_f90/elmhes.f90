	SUBROUTINE elmhes(a)
	USE nrtype; USE nrutil, ONLY : assert_eq,imaxloc,outerprod,swap
	IMPLICIT NONE
	REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
	INTEGER(I4B) :: i,m,n
	REAL(SP) :: x
	REAL(SP), DIMENSION(size(a,1)) :: y
	n=assert_eq(size(a,1),size(a,2),'elmhes')
	do m=2,n-1
		i=imaxloc(abs(a(m:n,m-1)))+m-1
		x=a(i,m-1)
		if (i /= m) then
			call swap(a(i,m-1:n),a(m,m-1:n))
			call swap(a(:,i),a(:,m))
		end if
		if (x /= 0.0) then
			y(m+1:n)=a(m+1:n,m-1)/x
			a(m+1:n,m-1)=y(m+1:n)
			a(m+1:n,m:n)=a(m+1:n,m:n)-outerprod(y(m+1:n),a(m,m:n))
			a(:,m)=a(:,m)+matmul(a(:,m+1:n),y(m+1:n))
		end if
	end do
	END SUBROUTINE elmhes
