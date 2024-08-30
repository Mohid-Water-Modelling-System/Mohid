	SUBROUTINE mpmul(w,u,v,n,m)
	USE nrtype; USE nrutil, ONLY : nrerror
	USE nr, ONLY : realft
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: n,m
	CHARACTER(1), DIMENSION(:), INTENT(IN) :: u,v
	CHARACTER(1), DIMENSION(:), INTENT(OUT) :: w
!	The logical dimensions are: CHARACTER(1) :: w(n+m),u(n),v(m)
	REAL(DP), PARAMETER :: RX=256.0
	INTEGER(I4B) :: j,mn,nn
	REAL(DP) :: cy,t
	REAL(DP), DIMENSION(:), ALLOCATABLE :: a,b,tb
	mn=max(m,n)
	nn=1
	do
		if (nn >= mn) exit
		nn=nn+nn
	end do
	nn=nn+nn
	allocate(a(nn),b(nn),tb((nn-1)/2))
	a(1:n)=ichar(u(1:n))
	a(n+1:nn)=0.0
	b(1:m)=ichar(v(1:m))
	b(m+1:nn)=0.0
	call realft(a(1:nn),1)
	call realft(b(1:nn),1)
	b(1)=b(1)*a(1)
	b(2)=b(2)*a(2)
	tb=b(3:nn:2)
	b(3:nn:2)=tb*a(3:nn:2)-b(4:nn:2)*a(4:nn:2)
	b(4:nn:2)=tb*a(4:nn:2)+b(4:nn:2)*a(3:nn:2)
	call realft(b(1:nn),-1)
	b(:)=b(:)/(nn/2)
	cy=0.0
	do j=nn,1,-1
		t=b(j)+cy+0.5_dp
		b(j)=mod(t,RX)
		cy=int(t/RX)
	end do
	if (cy >= RX) call nrerror('mpmul: sanity check failed in fftmul')
	w(1)=char(int(cy))
	w(2:(n+m))=char(int(b(1:(n+m-1))))
	deallocate(a,b,tb)
	END SUBROUTINE mpmul
