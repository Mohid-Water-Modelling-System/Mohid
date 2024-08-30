	PROGRAM fredex
	USE nrtype; USE nrutil, ONLY : arth
	USE nr, ONLY : quadmx,ludcmp,lubksb
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: N=40
	INTEGER(I4B) :: j
	INTEGER(I4B), DIMENSION(N) :: indx
	REAL(SP) :: d
	REAL(SP), DIMENSION(N) :: g,x
	REAL(SP), DIMENSION(N,N) :: a
	call quadmx(a)
	call ludcmp(a,indx,d)
	x(:)=arth(0,1,n)*PI/(n-1)
	g(:)=sin(x(:))
	call lubksb(a,indx,g)
	do j=1,n
		write (*,*) j,x(j),g(j)
	end do
	write (*,*) 'normal completion'
	END PROGRAM fredex
