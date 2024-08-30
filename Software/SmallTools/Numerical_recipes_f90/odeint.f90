MODULE ode_path
	USE nrtype
	INTEGER(I4B) :: nok,nbad,kount
	LOGICAL(LGT), SAVE :: save_steps=.false.
	REAL(SP) :: dxsav
	REAL(SP), DIMENSION(:), POINTER :: xp
	REAL(SP), DIMENSION(:,:), POINTER :: yp
END MODULE ode_path

	SUBROUTINE odeint(ystart,x1,x2,eps,h1,hmin,derivs,rkqs)
	USE nrtype; USE nrutil, ONLY : nrerror,reallocate
	USE ode_path
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: ystart
	REAL(SP), INTENT(IN) :: x1,x2,eps,h1,hmin
	INTERFACE
		SUBROUTINE derivs(x,y,dydx)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP), DIMENSION(:), INTENT(IN) :: y
		REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
		END SUBROUTINE derivs
!BL
		SUBROUTINE rkqs(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
		USE nrtype
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
		END SUBROUTINE rkqs
	END INTERFACE
	REAL(SP), PARAMETER :: TINY=1.0e-30_sp
	INTEGER(I4B), PARAMETER :: MAXSTP=10000
	INTEGER(I4B) :: nstp
	REAL(SP) :: h,hdid,hnext,x,xsav
	REAL(SP), DIMENSION(size(ystart)) :: dydx,y,yscal
	x=x1
	h=sign(h1,x2-x1)
	nok=0
	nbad=0
	kount=0
	y(:)=ystart(:)
	if (save_steps) then
		xsav=x-2.0_sp*dxsav
		nullify(xp,yp)
		allocate(xp(256))
		allocate(yp(size(ystart),size(xp)))
	end if
	do nstp=1,MAXSTP
		call derivs(x,y,dydx)
		yscal(:)=abs(y(:))+abs(h*dydx(:))+TINY
		if (save_steps .and. (abs(x-xsav) > abs(dxsav))) &
			call save_a_step
		if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x
		call rkqs(y,dydx,x,h,eps,yscal,hdid,hnext,derivs)
		if (hdid == h) then
			nok=nok+1
		else
			nbad=nbad+1
		end if
		if ((x-x2)*(x2-x1) >= 0.0) then
			ystart(:)=y(:)
			if (save_steps) call save_a_step
			RETURN
		end if
		if (abs(hnext) < hmin)&
			call nrerror('stepsize smaller than minimum in odeint')
		h=hnext
	end do
	call nrerror('too many steps in odeint')
	CONTAINS
!BL
	SUBROUTINE save_a_step
	kount=kount+1
	if (kount > size(xp)) then
		xp=>reallocate(xp,2*size(xp))
		yp=>reallocate(yp,size(yp,1),size(xp))
	end if
	xp(kount)=x
	yp(:,kount)=y(:)
	xsav=x
	END SUBROUTINE save_a_step
	END SUBROUTINE odeint
