	FUNCTION rtsec(func,x1,x2,xacc)
	USE nrtype; USE nrutil, ONLY : nrerror,swap
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x1,x2,xacc
	REAL(SP) :: rtsec
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
	REAL(SP) :: dx,f,fl,xl
	fl=func(x1)
	f=func(x2)
	if (abs(fl) < abs(f)) then
		rtsec=x1
		xl=x2
		call swap(fl,f)
	else
		xl=x1
		rtsec=x2
	end if
	do j=1,MAXIT
		dx=(xl-rtsec)*f/(f-fl)
		xl=rtsec
		fl=f
		rtsec=rtsec+dx
		f=func(rtsec)
		if (abs(dx) < xacc .or. f == 0.0) RETURN
	end do
	call nrerror('rtsec: exceed maximum iterations')
	END FUNCTION rtsec
