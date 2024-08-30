	FUNCTION zriddr(func,x1,x2,xacc)
	USE nrtype; USE nrutil, ONLY : nrerror
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x1,x2,xacc
	REAL(SP) :: zriddr
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER(I4B), PARAMETER :: MAXIT=60
	REAL(SP), PARAMETER :: UNUSED=-1.11e30_sp
	INTEGER(I4B) :: j
	REAL(SP) :: fh,fl,fm,fnew,s,xh,xl,xm,xnew
	fl=func(x1)
	fh=func(x2)
	if ((fl > 0.0 .and. fh < 0.0) .or. (fl < 0.0 .and. fh > 0.0)) then
		xl=x1
		xh=x2
		zriddr=UNUSED
		do j=1,MAXIT
			xm=0.5_sp*(xl+xh)
			fm=func(xm)
			s=sqrt(fm**2-fl*fh)
			if (s == 0.0) RETURN
			xnew=xm+(xm-xl)*(sign(1.0_sp,fl-fh)*fm/s)
			if (abs(xnew-zriddr) <= xacc) RETURN
			zriddr=xnew
			fnew=func(zriddr)
			if (fnew == 0.0) RETURN
			if (sign(fm,fnew) /= fm) then
				xl=xm
				fl=fm
				xh=zriddr
				fh=fnew
			else if (sign(fl,fnew) /= fl) then
				xh=zriddr
				fh=fnew
			else if (sign(fh,fnew) /= fh) then
				xl=zriddr
				fl=fnew
			else
				call nrerror('zriddr: never get here')
			end if
			if (abs(xh-xl) <= xacc) RETURN
		end do
		call nrerror('zriddr: exceeded maximum iterations')
	else if (fl == 0.0) then
		zriddr=x1
	else if (fh == 0.0) then
		zriddr=x2
	else
		call nrerror('zriddr: root must be bracketed')
	end if
	END FUNCTION zriddr
