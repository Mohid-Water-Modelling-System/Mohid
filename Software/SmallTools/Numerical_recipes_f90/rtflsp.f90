	FUNCTION rtflsp(func,x1,x2,xacc)
	USE nrtype; USE nrutil, ONLY : nrerror,swap
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x1,x2,xacc
	REAL(SP) :: rtflsp
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER(I4B), PARAMETER :: MAXIT=30
	INTEGER(I4B) :: j
	REAL(SP) :: del,dx,f,fh,fl,xh,xl
	fl=func(x1)
	fh=func(x2)
	if ((fl > 0.0 .and. fh > 0.0) .or. &
		(fl < 0.0 .and. fh < 0.0)) call &
		nrerror('rtflsp: root must be bracketed between arguments')
	if (fl < 0.0) then
		xl=x1
		xh=x2
	else
		xl=x2
		xh=x1
		call swap(fl,fh)
	end if
	dx=xh-xl
	do j=1,MAXIT
		rtflsp=xl+dx*fl/(fl-fh)
		f=func(rtflsp)
		if (f < 0.0) then
			del=xl-rtflsp
			xl=rtflsp
			fl=f
		else
			del=xh-rtflsp
			xh=rtflsp
			fh=f
		end if
		dx=xh-xl
		if (abs(del) < xacc .or. f == 0.0) RETURN
	end do
	call nrerror('rtflsp exceed maximum iterations')
	END FUNCTION rtflsp
