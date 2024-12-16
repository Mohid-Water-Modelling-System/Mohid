	FUNCTION golden(ax,bx,cx,func,tol,xmin)
	USE nrtype
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: ax,bx,cx,tol
	REAL(SP), INTENT(OUT) :: xmin
	REAL(SP) :: golden
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: func
		END FUNCTION func
	END INTERFACE
	REAL(SP), PARAMETER :: R=0.61803399_sp,C=1.0_sp-R
	REAL(SP) :: f1,f2,x0,x1,x2,x3
	x0=ax
	x3=cx
	if (abs(cx-bx) > abs(bx-ax)) then
		x1=bx
		x2=bx+C*(cx-bx)
	else
		x2=bx
		x1=bx-C*(bx-ax)
	end if
	f1=func(x1)
	f2=func(x2)
	do
		if (abs(x3-x0) <= tol*(abs(x1)+abs(x2))) exit
		if (f2 < f1) then
			call shft3(x0,x1,x2,R*x2+C*x3)
			call shft2(f1,f2,func(x2))
		else
			call shft3(x3,x2,x1,R*x1+C*x0)
			call shft2(f2,f1,func(x1))
		end if
	end do
	if (f1 < f2) then
		golden=f1
		xmin=x1
	else
		golden=f2
		xmin=x2
	end if
	CONTAINS
!BL
	SUBROUTINE shft2(a,b,c)
	REAL(SP), INTENT(OUT) :: a
	REAL(SP), INTENT(INOUT) :: b
	REAL(SP), INTENT(IN) :: c
	a=b
	b=c
	END SUBROUTINE shft2
!BL
	SUBROUTINE shft3(a,b,c,d)
	REAL(SP), INTENT(OUT) :: a
	REAL(SP), INTENT(INOUT) :: b,c
	REAL(SP), INTENT(IN) :: d
	a=b
	b=c
	c=d
	END SUBROUTINE shft3
	END FUNCTION golden
