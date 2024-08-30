	FUNCTION dfridr(func,x,h,err)
	USE nrtype; USE nrutil, ONLY : assert,geop,iminloc
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x,h
	REAL(SP), INTENT(OUT) :: err
	REAL(SP) :: dfridr
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER(I4B),PARAMETER :: NTAB=10
	REAL(SP), PARAMETER :: CON=1.4_sp,CON2=CON*CON,BIG=huge(x),SAFE=2.0
	INTEGER(I4B) :: ierrmin,i,j
	REAL(SP) :: hh
	REAL(SP), DIMENSION(NTAB-1) :: errt,fac
	REAL(SP), DIMENSION(NTAB,NTAB) :: a
	call assert(h /= 0.0, 'dfridr arg')
	hh=h
	a(1,1)=(func(x+hh)-func(x-hh))/(2.0_sp*hh)
	err=BIG
	fac(1:NTAB-1)=geop(CON2,CON2,NTAB-1)
	do i=2,NTAB
		hh=hh/CON
		a(1,i)=(func(x+hh)-func(x-hh))/(2.0_sp*hh)
		do j=2,i
			a(j,i)=(a(j-1,i)*fac(j-1)-a(j-1,i-1))/(fac(j-1)-1.0_sp)
		end do
		errt(1:i-1)=max(abs(a(2:i,i)-a(1:i-1,i)),abs(a(2:i,i)-a(1:i-1,i-1)))
		ierrmin=iminloc(errt(1:i-1))
		if (errt(ierrmin) <= err) then
			err=errt(ierrmin)
			dfridr=a(1+ierrmin,i)
		end if
		if (abs(a(i,i)-a(i-1,i-1)) >= SAFE*err) RETURN
	end do
	END FUNCTION dfridr
