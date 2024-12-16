	SUBROUTINE rzextr(iest,xest,yest,yz,dy)
	USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: iest
	REAL(SP), INTENT(IN) :: xest
	REAL(SP), DIMENSION(:), INTENT(IN) :: yest
	REAL(SP), DIMENSION(:), INTENT(OUT) :: yz,dy
	INTEGER(I4B), PARAMETER :: IEST_MAX=16
	INTEGER(I4B) :: k,nv
	INTEGER(I4B), SAVE :: nvold=-1
	REAL(SP), DIMENSION(size(yz)) :: yy,v,c,b,b1,ddy
	REAL(SP), DIMENSION(:,:), ALLOCATABLE, SAVE :: d
	REAL(SP), DIMENSION(IEST_MAX), SAVE :: fx,x
	nv=assert_eq(size(yz),size(dy),size(yest),'rzextr')
	if (iest > IEST_MAX) call &
		nrerror('rzextr: probable misuse, too much extrapolation')
	if (nv /= nvold) then
		if (allocated(d)) deallocate(d)
		allocate(d(nv,IEST_MAX))
		nvold=nv
	end if
	x(iest)=xest
	if (iest == 1) then
		yz=yest
		d(:,1)=yest
		dy=yest
	else
		fx(2:iest)=x(iest-1:1:-1)/xest
		yy=yest
		v=d(1:nv,1)
		c=yy
		d(1:nv,1)=yy
		do k=2,iest
			b1=fx(k)*v
			b=b1-c
			where (b /= 0.0)
				b=(c-v)/b
				ddy=c*b
				c=b1*b
			elsewhere
				ddy=v
			end where
			if (k /= iest) v=d(1:nv,k)
			d(1:nv,k)=ddy
			yy=yy+ddy
		end do
		dy=ddy
		yz=yy
	end if
	END SUBROUTINE rzextr
