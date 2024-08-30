	SUBROUTINE pzextr(iest,xest,yest,yz,dy)
	USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: iest
	REAL(SP), INTENT(IN) :: xest
	REAL(SP), DIMENSION(:), INTENT(IN) :: yest
	REAL(SP), DIMENSION(:), INTENT(OUT) :: yz,dy
	INTEGER(I4B), PARAMETER :: IEST_MAX=16
	INTEGER(I4B) :: j,nv
	INTEGER(I4B), SAVE :: nvold=-1
	REAL(SP) :: delta,f1,f2
	REAL(SP), DIMENSION(size(yz)) :: d,tmp,q
	REAL(SP), DIMENSION(IEST_MAX), SAVE :: x
	REAL(SP), DIMENSION(:,:), ALLOCATABLE, SAVE :: qcol
	nv=assert_eq(size(yz),size(yest),size(dy),'pzextr')
	if (iest > IEST_MAX) call &
		nrerror('pzextr: probable misuse, too much extrapolation')
	if (nv /= nvold) then
		if (allocated(qcol)) deallocate(qcol)
		allocate(qcol(nv,IEST_MAX))
		nvold=nv
	end if
	x(iest)=xest
	dy(:)=yest(:)
	yz(:)=yest(:)
	if (iest == 1) then
		qcol(:,1)=yest(:)
	else
		d(:)=yest(:)
		do j=1,iest-1
			delta=1.0_sp/(x(iest-j)-xest)
			f1=xest*delta
			f2=x(iest-j)*delta
			q(:)=qcol(:,j)
			qcol(:,j)=dy(:)
			tmp(:)=d(:)-q(:)
			dy(:)=f1*tmp(:)
			d(:)=f2*tmp(:)
			yz(:)=yz(:)+dy(:)
		end do
		qcol(:,iest)=dy(:)
	end if
	END SUBROUTINE pzextr
