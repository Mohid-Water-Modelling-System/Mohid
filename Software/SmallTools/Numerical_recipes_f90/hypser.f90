	SUBROUTINE hypser(a,b,c,z,series,deriv)
	USE nrtype; USE nrutil, ONLY : nrerror
	IMPLICIT NONE
	COMPLEX(SPC), INTENT(IN) :: a,b,c,z
	COMPLEX(SPC), INTENT(OUT) :: series,deriv
	INTEGER(I4B) :: n
	INTEGER(I4B), PARAMETER :: MAXIT=1000
	COMPLEX(SPC) :: aa,bb,cc,fac,temp
	deriv=cmplx(0.0_sp,0.0_sp,kind=spc)
	fac=cmplx(1.0_sp,0.0_sp,kind=spc)
	temp=fac
	aa=a
	bb=b
	cc=c
	do n=1,MAXIT
		fac=((aa*bb)/cc)*fac
		deriv=deriv+fac
		fac=fac*z/n
		series=temp+fac
		if (series == temp) RETURN
		temp=series
		aa=aa+1.0
		bb=bb+1.0
		cc=cc+1.0
	end do
	call nrerror('hypser: convergence failure')
	END SUBROUTINE hypser
