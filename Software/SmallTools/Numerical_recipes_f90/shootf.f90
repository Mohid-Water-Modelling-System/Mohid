!	FUNCTION shootf(v) is named "funcv" for use with "newt"
	FUNCTION funcv(v)
	USE nrtype
	USE nr, ONLY : odeint,rkqs
	USE sphfpt_caller, ONLY : x1,x2,xf,nn2; USE ode_path, ONLY : xp,yp
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: v
	REAL(SP), DIMENSION(size(v)) :: funcv
	REAL(SP), PARAMETER :: EPS=1.0e-6_sp
	REAL(SP) :: h1,hmin
	REAL(SP), DIMENSION(size(v)) :: f1,f2,y
	INTERFACE
		SUBROUTINE derivs(x,y,dydx)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP), DIMENSION(:), INTENT(IN) :: y
		REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
		END SUBROUTINE derivs
!BL
		SUBROUTINE load1(x1,v1,y)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x1
		REAL(SP), DIMENSION(:), INTENT(IN) :: v1
		REAL(SP), DIMENSION(:), INTENT(OUT) :: y
		END SUBROUTINE load1
!BL
		SUBROUTINE load2(x2,v2,y)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x2
		REAL(SP), DIMENSION(:), INTENT(IN) :: v2
		REAL(SP), DIMENSION(:), INTENT(OUT) :: y
		END SUBROUTINE load2
!BL
		SUBROUTINE score(x2,y,f)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x2
		REAL(SP), DIMENSION(:), INTENT(IN) :: y
		REAL(SP), DIMENSION(:), INTENT(OUT) :: f
		END SUBROUTINE score
	END INTERFACE
	h1=(x2-x1)/100.0_sp
	hmin=0.0
	call load1(x1,v,y)
	if (associated(xp)) deallocate(xp,yp)
	call odeint(y,x1,xf,EPS,h1,hmin,derivs,rkqs)
	call score(xf,y,f1)
	call load2(x2,v(nn2+1:),y)
	call odeint(y,x2,xf,EPS,h1,hmin,derivs,rkqs)
	call score(xf,y,f2)
	funcv(:)=f1(:)-f2(:)
	END FUNCTION funcv
