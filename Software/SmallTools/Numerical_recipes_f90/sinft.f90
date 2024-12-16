	SUBROUTINE sinft(y)
	USE nrtype; USE nrutil, ONLY : assert,cumsum,zroots_unity
	USE nr, ONLY : realft
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: y
	REAL(SP), DIMENSION(size(y)/2+1) :: wi
	REAL(SP), DIMENSION(size(y)/2) :: y1,y2
	INTEGER(I4B) :: n,nh
	n=size(y)
	call assert(iand(n,n-1)==0, 'n must be a power of 2 in sinft')
	nh=n/2
	wi=aimag(zroots_unity(n+n,nh+1))
	y(1)=0.0
	y1=wi(2:nh+1)*(y(2:nh+1)+y(n:nh+1:-1))
	y2=0.5_sp*(y(2:nh+1)-y(n:nh+1:-1))
	y(2:nh+1)=y1+y2
	y(n:nh+1:-1)=y1-y2
	call realft(y,+1)
	y(1)=0.5_sp*y(1)
	y(2)=0.0
	y1=cumsum(y(1:n-1:2))
	y(1:n-1:2)=y(2:n:2)
	y(2:n:2)=y1
	END SUBROUTINE sinft
