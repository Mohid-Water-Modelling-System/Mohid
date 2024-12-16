	SUBROUTINE zrhqr(a,rtr,rti)
	USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
	USE nr, ONLY : balanc,hqr,indexx
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: a
	REAL(SP), DIMENSION(:), INTENT(OUT) :: rtr,rti
	INTEGER(I4B) :: k,m
	INTEGER(I4B), DIMENSION(size(rtr)) :: indx
	REAL(SP), DIMENSION(size(a)-1,size(a)-1) :: hess
	m=assert_eq(size(rtr),size(rti),size(a)-1,'zrhqr')
	if (a(m+1) == 0.0) call &
		nrerror('zrhqr: Last value of array a must not be 0')
	hess(1,:)=-a(m:1:-1)/a(m+1)
	hess(2:m,:)=0.0
	do k=1,m-1
		hess(k+1,k)=1.0
	end do
	call balanc(hess)
	call hqr(hess,rtr,rti)
	call indexx(rtr,indx)
	rtr=rtr(indx)
	rti=rti(indx)
	END SUBROUTINE zrhqr
