	SUBROUTINE pade(cof,resid)
	USE nrtype
	USE nr, ONLY : lubksb,ludcmp,mprove
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), INTENT(INOUT) :: cof
	REAL(SP), INTENT(OUT) :: resid
	INTEGER(I4B) :: k,n
	INTEGER(I4B), DIMENSION((size(cof)-1)/2) :: indx
	REAL(SP), PARAMETER :: BIG=1.0e30_sp
	REAL(SP) :: d,rr,rrold
	REAL(SP), DIMENSION((size(cof)-1)/2) :: x,y,z
	REAL(SP), DIMENSION((size(cof)-1)/2,(size(cof)-1)/2) :: q,qlu
	n=(size(cof)-1)/2
	x=cof(n+2:2*n+1)
	y=x
	do k=1,n
		q(:,k)=cof(n+2-k:2*n+1-k)
	end do
	qlu=q
	call ludcmp(qlu,indx,d)
	call lubksb(qlu,indx,x)
	rr=BIG
	do
		rrold=rr
		z=x
		call mprove(q,qlu,indx,y,x)
		rr=sum((z-x)**2)
		if (rr >= rrold) exit
	end do
	resid=sqrt(rrold)
	do k=1,n
		y(k)=cof(k+1)-dot_product(z(1:k),cof(k:1:-1))
	end do
	cof(2:n+1)=y
	cof(n+2:2*n+1)=-z
	END SUBROUTINE pade
