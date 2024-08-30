	SUBROUTINE cosft2(y,isign)
	USE nrtype; USE nrutil, ONLY : assert,cumsum,zroots_unity
	USE nr, ONLY : realft
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: y
	INTEGER(I4B), INTENT(IN) :: isign
	COMPLEX(SPC), DIMENSION(size(y)) :: w
	REAL(SP), DIMENSION(size(y)/2) :: y1,y2
	REAL(SP) :: ytemp
	INTEGER(I4B) :: n,nh
	n=size(y)
	call assert(iand(n,n-1)==0, 'n must be a power of 2 in cosft2')
	nh=n/2
	w=zroots_unity(4*n,n)
	if (isign == 1) then
		y1=0.5_sp*(y(1:nh)+y(n:nh+1:-1))
		y2=aimag(w(2:n:2))*(y(1:nh)-y(n:nh+1:-1))
		y(1:nh)=y1+y2
		y(n:nh+1:-1)=y1-y2
		call realft(y,1)
		y1(1:nh-1)=y(3:n-1:2)*real(w(3:n-1:2)) &
			-y(4:n:2)*aimag(w(3:n-1:2))
		y2(1:nh-1)=y(4:n:2)*real(w(3:n-1:2)) &
			+y(3:n-1:2)*aimag(w(3:n-1:2))
		y(3:n-1:2)=y1(1:nh-1)
		y(4:n:2)=y2(1:nh-1)
		ytemp=0.5_sp*y(2)
		y(n-2:2:-2)=cumsum(y(n:4:-2),ytemp)
		y(n)=ytemp
	else if (isign == -1) then
		ytemp=y(n)
		y(4:n:2)=y(2:n-2:2)-y(4:n:2)
		y(2)=2.0_sp*ytemp
		y1(1:nh-1)=y(3:n-1:2)*real(w(3:n-1:2)) &
			+y(4:n:2)*aimag(w(3:n-1:2))
		y2(1:nh-1)=y(4:n:2)*real(w(3:n-1:2)) &
			-y(3:n-1:2)*aimag(w(3:n-1:2))
		y(3:n-1:2)=y1(1:nh-1)
		y(4:n:2)=y2(1:nh-1)
		call realft(y,-1)
		y1=y(1:nh)+y(n:nh+1:-1)
		y2=(0.5_sp/aimag(w(2:n:2)))*(y(1:nh)-y(n:nh+1:-1))
		y(1:nh)=0.5_sp*(y1+y2)
		y(n:nh+1:-1)=0.5_sp*(y1-y2)
	end if
	END SUBROUTINE cosft2
