	SUBROUTINE rkqs(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
	USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
	USE nr, ONLY : rkck
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: y
	REAL(SP), DIMENSION(:), INTENT(IN) :: dydx,yscal
	REAL(SP), INTENT(INOUT) :: x
	REAL(SP), INTENT(IN) :: htry,eps
	REAL(SP), INTENT(OUT) :: hdid,hnext
	INTERFACE
		SUBROUTINE derivs(x,y,dydx)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP), DIMENSION(:), INTENT(IN) :: y
		REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
		END SUBROUTINE derivs
	END INTERFACE
	INTEGER(I4B) :: ndum
	REAL(SP) :: errmax,h,htemp,xnew
	REAL(SP), DIMENSION(size(y)) :: yerr,ytemp
	REAL(SP), PARAMETER :: SAFETY=0.9_sp,PGROW=-0.2_sp,PSHRNK=-0.25_sp,&
		ERRCON=1.89e-4
	ndum=assert_eq(size(y),size(dydx),size(yscal),'rkqs')
	h=htry
	do
		call rkck(y,dydx,x,h,ytemp,yerr,derivs)
		errmax=maxval(abs(yerr(:)/yscal(:)))/eps
		if (errmax <= 1.0) exit
		htemp=SAFETY*h*(errmax**PSHRNK)
		h=sign(max(abs(htemp),0.1_sp*abs(h)),h)
		xnew=x+h
		if (xnew == x) call nrerror('stepsize underflow in rkqs')
	end do
	if (errmax > ERRCON) then
		hnext=SAFETY*h*(errmax**PGROW)
	else
		hnext=5.0_sp*h
	end if
	hdid=h
	x=x+h
	y(:)=ytemp(:)
	END SUBROUTINE rkqs
