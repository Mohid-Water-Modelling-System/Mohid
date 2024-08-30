	SUBROUTINE rkck(y,dydx,x,h,yout,yerr,derivs)
	USE nrtype; USE nrutil, ONLY : assert_eq
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: y,dydx
	REAL(SP), INTENT(IN) :: x,h
	REAL(SP), DIMENSION(:), INTENT(OUT) :: yout,yerr
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
	REAL(SP), DIMENSION(size(y)) :: ak2,ak3,ak4,ak5,ak6,ytemp
	REAL(SP), PARAMETER :: A2=0.2_sp,A3=0.3_sp,A4=0.6_sp,A5=1.0_sp,&
		A6=0.875_sp,B21=0.2_sp,B31=3.0_sp/40.0_sp,B32=9.0_sp/40.0_sp,&
		B41=0.3_sp,B42=-0.9_sp,B43=1.2_sp,B51=-11.0_sp/54.0_sp,&
		B52=2.5_sp,B53=-70.0_sp/27.0_sp,B54=35.0_sp/27.0_sp,&
		B61=1631.0_sp/55296.0_sp,B62=175.0_sp/512.0_sp,&
		B63=575.0_sp/13824.0_sp,B64=44275.0_sp/110592.0_sp,&
		B65=253.0_sp/4096.0_sp,C1=37.0_sp/378.0_sp,&
		C3=250.0_sp/621.0_sp,C4=125.0_sp/594.0_sp,&
		C6=512.0_sp/1771.0_sp,DC1=C1-2825.0_sp/27648.0_sp,&
		DC3=C3-18575.0_sp/48384.0_sp,DC4=C4-13525.0_sp/55296.0_sp,&
		DC5=-277.0_sp/14336.0_sp,DC6=C6-0.25_sp
	ndum=assert_eq(size(y),size(dydx),size(yout),size(yerr),'rkck')
	ytemp=y+B21*h*dydx
	call derivs(x+A2*h,ytemp,ak2)
	ytemp=y+h*(B31*dydx+B32*ak2)
	call derivs(x+A3*h,ytemp,ak3)
	ytemp=y+h*(B41*dydx+B42*ak2+B43*ak3)
	call derivs(x+A4*h,ytemp,ak4)
	ytemp=y+h*(B51*dydx+B52*ak2+B53*ak3+B54*ak4)
	call derivs(x+A5*h,ytemp,ak5)
	ytemp=y+h*(B61*dydx+B62*ak2+B63*ak3+B64*ak4+B65*ak5)
	call derivs(x+A6*h,ytemp,ak6)
	yout=y+h*(C1*dydx+C3*ak3+C4*ak4+C6*ak6)
	yerr=h*(DC1*dydx+DC3*ak3+DC4*ak4+DC5*ak5+DC6*ak6)
	END SUBROUTINE rkck
