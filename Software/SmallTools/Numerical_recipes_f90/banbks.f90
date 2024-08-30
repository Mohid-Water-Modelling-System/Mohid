	SUBROUTINE banbks(a,m1,m2,al,indx,b)
	USE nrtype; USE nrutil, ONLY : assert_eq,swap
	IMPLICIT NONE
	REAL(SP), DIMENSION(:,:), INTENT(IN) :: a,al
	INTEGER(I4B), INTENT(IN) :: m1,m2
	INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: b
	INTEGER(I4B) :: i,k,l,mdum,mm,n
	n=assert_eq(size(a,1),size(al,1),size(b),size(indx),'banbks: n')
	mm=assert_eq(size(a,2),m1+m2+1,'banbks: mm')
	mdum=assert_eq(size(al,2),m1,'banbks: mdum')
	do k=1,n
		l=min(n,m1+k)
		i=indx(k)
		if (i /= k) call swap(b(i),b(k))
		b(k+1:l)=b(k+1:l)-al(k,1:l-k)*b(k)
	end do
	do i=n,1,-1
		l=min(mm,n-i+1)
		b(i)=(b(i)-dot_product(a(i,2:l),b(1+i:i+l-1)))/a(i,1)
	end do
	END SUBROUTINE banbks
