	SUBROUTINE zbrac(func,x1,x2,succes)
	USE nrtype; USE nrutil, ONLY : nrerror
	IMPLICIT NONE
	REAL(SP), INTENT(INOUT) :: x1,x2
	LOGICAL(LGT), INTENT(OUT) :: succes
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER(I4B), PARAMETER :: NTRY=50
	REAL(SP), PARAMETER :: FACTOR=1.6_sp
	INTEGER(I4B) :: j
	REAL(SP) :: f1,f2
	if (x1 == x2) call nrerror('zbrac: you have to guess an initial range')
	f1=func(x1)
	f2=func(x2)
	succes=.true.
	do j=1,NTRY
		if ((f1 > 0.0 .and. f2 < 0.0) .or. &
			(f1 < 0.0 .and. f2 > 0.0)) RETURN
		if (abs(f1) < abs(f2)) then
			x1=x1+FACTOR*(x1-x2)
			f1=func(x1)
		else
			x2=x2+FACTOR*(x2-x1)
			f2=func(x2)
		end if
	end do
	succes=.false.
	END SUBROUTINE zbrac
