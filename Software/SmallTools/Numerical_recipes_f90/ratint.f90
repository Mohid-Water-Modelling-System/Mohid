	SUBROUTINE ratint(xa,ya,x,y,dy)
	USE nrtype; USE nrutil, ONLY : assert_eq,iminloc,nrerror
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: xa,ya
	REAL(SP), INTENT(IN) :: x
	REAL(SP), INTENT(OUT) :: y,dy
	INTEGER(I4B) :: m,n,ns
	REAL(SP), DIMENSION(size(xa)) :: c,d,dd,h,t
	REAL(SP), PARAMETER :: TINY=1.0e-25_sp
	n=assert_eq(size(xa),size(ya),'ratint')
	h=xa-x
	ns=iminloc(abs(h))
	y=ya(ns)
	if (x == xa(ns)) then
		dy=0.0
		RETURN
	end if
	c=ya
	d=ya+TINY
	ns=ns-1
	do m=1,n-1
		t(1:n-m)=(xa(1:n-m)-x)*d(1:n-m)/h(1+m:n)
		dd(1:n-m)=t(1:n-m)-c(2:n-m+1)
		if (any(dd(1:n-m) == 0.0)) &
			call nrerror('failure in ratint')
		dd(1:n-m)=(c(2:n-m+1)-d(1:n-m))/dd(1:n-m)
		d(1:n-m)=c(2:n-m+1)*dd(1:n-m)
		c(1:n-m)=t(1:n-m)*dd(1:n-m)
		if (2*ns < n-m) then
			dy=c(ns+1)
		else
			dy=d(ns)
			ns=ns-1
		end if
		y=y+dy
	end do
	END SUBROUTINE ratint
