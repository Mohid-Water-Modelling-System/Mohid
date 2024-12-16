	SUBROUTINE stoerm(y,d2y,xs,htot,nstep,yout,derivs)
	USE nrtype; USE nrutil, ONLY : assert_eq
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: y,d2y
	REAL(SP), INTENT(IN) :: xs,htot
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
	INTEGER(I4B) :: neqn,neqn1,nn,nv
	REAL(SP) :: h,h2,halfh,x
	REAL(SP), DIMENSION(size(y)) :: ytemp
	nv=assert_eq(size(y),size(d2y),size(yout),'stoerm')
	neqn=nv/2
	neqn1=neqn+1
	h=htot/nstep
	halfh=0.5_sp*h
	ytemp(neqn1:nv)=h*(y(neqn1:nv)+halfh*d2y(1:neqn))
	ytemp(1:neqn)=y(1:neqn)+ytemp(neqn1:nv)
	x=xs+h
	call derivs(x,ytemp,yout)
	h2=h*h
	do nn=2,nstep
		ytemp(neqn1:nv)=ytemp(neqn1:nv)+h2*yout(1:neqn)
		ytemp(1:neqn)=ytemp(1:neqn)+ytemp(neqn1:nv)
		x=x+h
		call derivs(x,ytemp,yout)
	end do
	yout(neqn1:nv)=ytemp(neqn1:nv)/h+halfh*yout(1:neqn)
	yout(1:neqn)=ytemp(1:neqn)
	END SUBROUTINE stoerm
