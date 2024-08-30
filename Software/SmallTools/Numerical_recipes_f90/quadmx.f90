	SUBROUTINE quadmx(a)
	USE nrtype; USE nrutil, ONLY : arth,assert_eq,diagadd,outerprod
	USE nr, ONLY : wwghts,kermom
	USE kermom_info
	IMPLICIT NONE
	REAL(SP), DIMENSION(:,:), INTENT(OUT) :: a
	INTEGER(I4B) :: j,n
	REAL(SP) :: h,x
	REAL(SP), DIMENSION(size(a,1)) :: wt
	n=assert_eq(size(a,1),size(a,2),'quadmx')
	h=PI/(n-1)
	do j=1,n
		x=(j-1)*h
		kermom_x=x
		wt(:)=wwghts(n,h,kermom)
		a(j,:)=wt(:)
	end do
	wt(:)=cos(arth(0,1,n)*h)
	a(:,:)=a(:,:)*outerprod(wt(:),wt(:))
	call diagadd(a,1.0_sp)
	END SUBROUTINE quadmx
