	SUBROUTINE ddpoly(c,x,pd)
	USE nrtype; USE nrutil, ONLY : arth,cumprod,poly_term
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP), DIMENSION(:), INTENT(IN) :: c
	REAL(SP), DIMENSION(:), INTENT(OUT) :: pd
	INTEGER(I4B) :: i,nc,nd
	REAL(SP), DIMENSION(size(pd)) :: fac
	REAL(SP), DIMENSION(size(c)) :: d
	nc=size(c)
	nd=size(pd)
	d(nc:1:-1)=poly_term(c(nc:1:-1),x)
	do i=2,min(nd,nc)
		d(nc:i:-1)=poly_term(d(nc:i:-1),x)
	end do
	pd=d(1:nd)
	fac=cumprod(arth(1.0_sp,1.0_sp,nd))
	pd(3:nd)=fac(2:nd-1)*pd(3:nd)
	END SUBROUTINE ddpoly
