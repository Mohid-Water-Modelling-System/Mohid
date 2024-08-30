	SUBROUTINE cosft1(y)
	USE nrtype; USE nrutil, ONLY : assert,cumsum,zroots_unity
	USE nr, ONLY : realft
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: y
	COMPLEX(SPC), DIMENSION((size(y)-1)/2) :: w
	REAL(SP), DIMENSION((size(y)-1)/2-1) :: y1,y2
	REAL(SP) :: summ
	INTEGER(I4B) :: n,nh
	n=size(y)-1
	call assert(iand(n,n-1)==0, 'n must be a power of 2 in cosft1')
	nh=n/2
	w=zroots_unity(n+n,nh)
	summ=0.5_sp*(y(1)-y(n+1))
	y(1)=0.5_sp*(y(1)+y(n+1))
	y1=0.5_sp*(y(2:nh)+y(n:nh+2:-1))
	y2=y(2:nh)-y(n:nh+2:-1)
	summ=summ+sum(real(w(2:nh))*y2)
	y2=y2*aimag(w(2:nh))
	y(2:nh)=y1-y2
	y(n:nh+2:-1)=y1+y2
	call realft(y(1:n),1)
	y(n+1)=y(2)
	y(2)=summ
	y(2:n:2)=cumsum(y(2:n:2))
	END SUBROUTINE cosft1
