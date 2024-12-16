	FUNCTION expint(n,x)
	USE nrtype; USE nrutil, ONLY : arth,assert,nrerror
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: n
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: expint
	INTEGER(I4B), PARAMETER :: MAXIT=100
	REAL(SP), PARAMETER :: EPS=epsilon(x),BIG=huge(x)*EPS
	INTEGER(I4B) :: i,nm1
	REAL(SP) :: a,b,c,d,del,fact,h
	call assert(n >= 0, x >= 0.0, (x > 0.0 .or. n > 1), &
		'expint args')
	if (n == 0) then
		expint=exp(-x)/x
		RETURN
	end if
	nm1=n-1
	if (x == 0.0) then
		expint=1.0_sp/nm1
	else if (x > 1.0) then
		b=x+n
		c=BIG
		d=1.0_sp/b
		h=d
		do i=1,MAXIT
			a=-i*(nm1+i)
			b=b+2.0_sp
			d=1.0_sp/(a*d+b)
			c=b+a/c
			del=c*d
			h=h*del
			if (abs(del-1.0_sp) <= EPS) exit
		end do
		if (i > MAXIT) call nrerror('expint: continued fraction failed')
		expint=h*exp(-x)
	else
		if (nm1 /= 0) then
			expint=1.0_sp/nm1
		else
			expint=-log(x)-EULER
		end if
		fact=1.0
		do i=1,MAXIT
			fact=-fact*x/i
			if (i /= nm1) then
				del=-fact/(i-nm1)
			else
				del=fact*(-log(x)-EULER+sum(1.0_sp/arth(1,1,nm1)))
			end if
			expint=expint+del
			if (abs(del) < abs(expint)*EPS) exit
		end do
		if (i > MAXIT) call nrerror('expint: series failed')
	end if
	END FUNCTION expint
