	FUNCTION polcof(xa,ya)
	USE nrtype; USE nrutil, ONLY : assert_eq,iminloc
	USE nr, ONLY : polint
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: xa,ya
	REAL(SP), DIMENSION(size(xa)) :: polcof
	INTEGER(I4B) :: j,k,m,n
	REAL(SP) :: dy
	REAL(SP), DIMENSION(size(xa)) :: x,y
	n=assert_eq(size(xa),size(ya),'polcof')
	x=xa
	y=ya
	do j=1,n
		m=n+1-j
		call polint(x(1:m),y(1:m),0.0_sp,polcof(j),dy)
		k=iminloc(abs(x(1:m)))
		where (x(1:m) /= 0.0) y(1:m)=(y(1:m)-polcof(j))/x(1:m)
		y(k:m-1)=y(k+1:m)
		x(k:m-1)=x(k+1:m)
	end do
	END FUNCTION polcof
