	SUBROUTINE mppi(n)
	USE nrtype
	USE nr, ONLY : mp2dfr,mpinv,mpmul,mpsqrt
	USE mpops, ONLY : mpadd,mplsh,mpmov,mpsdv
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: n
	INTEGER(I4B), PARAMETER :: IAOFF=48
	INTEGER(I4B) :: ir,j,m
	CHARACTER(1), DIMENSION(n) :: sx,sxi
	CHARACTER(1), DIMENSION(2*n) :: t,y
	CHARACTER(1), DIMENSION(3*n) :: s
	CHARACTER(1), DIMENSION(n+1) :: x,bigpi
	t(1)=char(2)
	t(2:n)=char(0)
	call mpsqrt(x,x,t,n,n)
	call mpadd(bigpi,t,x,n)
	call mplsh(bigpi,n)
	call mpsqrt(sx,sxi,x,n,n)
	call mpmov(y,sx,n)
	do
		call mpadd(x,sx,sxi,n)
		call mpsdv(x,x(2:),n,2,ir)
		call mpsqrt(sx,sxi,x,n,n)
		call mpmul(t,y,sx,n,n)
		call mpadd(t(2:),t(2:),sxi,n)
		x(1)=char(ichar(x(1))+1)
		y(1)=char(ichar(y(1))+1)
		call mpinv(s,y,n,n)
		call mpmul(y,t(3:),s,n,n)
		call mplsh(y,n)
		call mpmul(t,x,s,n,n)
		m=mod(255+ichar(t(2)),256)
		if (abs(ichar(t(n+1))-m) > 1 .or. any(ichar(t(3:n)) /= m)) then
			call mpmul(s,bigpi,t(2:),n,n)
			call mpmov(bigpi,s(2:),n)
			cycle
		end if
		write (*,*) 'pi='
		s(1)=char(ichar(bigpi(1))+IAOFF)
		s(2)='.'
		call mp2dfr(bigpi(2:),s(3:),n-1,m)
		write (*,'(1x,64a1)') (s(j),j=1,m+1)
		RETURN
	end do
	END SUBROUTINE mppi
