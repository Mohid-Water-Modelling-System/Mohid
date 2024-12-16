	SUBROUTINE simpr(y,dydx,dfdx,dfdy,xs,htot,nstep,yout,derivs)
	USE nrtype; USE nrutil, ONLY : assert_eq,diagadd
	USE nr, ONLY : lubksb,ludcmp
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: xs,htot
	REAL(SP), DIMENSION(:), INTENT(IN) :: y,dydx,dfdx
	REAL(SP), DIMENSION(:,:), INTENT(IN) :: dfdy
	INTEGER(I4B), INTENT(IN) :: nstep
	REAL(SP), DIMENSION(:), INTENT(OUT) :: yout
	INTERFACE
		SUBROUTINE derivs(x,y,dydx)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP), DIMENSION(:), INTENT(IN) :: y
		REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
		END SUBROUTINE derivs
	END INTERFACE
	INTEGER(I4B) :: ndum,nn
	INTEGER(I4B), DIMENSION(size(y)) :: indx
	REAL(SP) :: d,h,x
	REAL(SP), DIMENSION(size(y)) :: del,ytemp
	REAL(SP), DIMENSION(size(y),size(y)) :: a
	ndum=assert_eq((/size(y),size(dydx),size(dfdx),size(dfdy,1),&
		size(dfdy,2),size(yout)/),'simpr')
	h=htot/nstep
	a(:,:)=-h*dfdy(:,:)
	call diagadd(a,1.0_sp)
	call ludcmp(a,indx,d)
	yout=h*(dydx+h*dfdx)
	call lubksb(a,indx,yout)
	del=yout
	ytemp=y+del
	x=xs+h
	call derivs(x,ytemp,yout)
	do nn=2,nstep
		yout=h*yout-del
		call lubksb(a,indx,yout)
		del=del+2.0_sp*yout
		ytemp=ytemp+del
		x=x+h
		call derivs(x,ytemp,yout)
	end do
	yout=h*yout-del
	call lubksb(a,indx,yout)
	yout=ytemp+yout
	END SUBROUTINE simpr
