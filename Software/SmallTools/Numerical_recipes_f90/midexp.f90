	SUBROUTINE midexp(funk,aa,bb,s,n)
	USE nrtype; USE nrutil, ONLY : arth
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: aa,bb
	REAL(SP), INTENT(INOUT) :: s
	INTEGER(I4B), INTENT(IN) :: n
	INTERFACE
		FUNCTION funk(x)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: funk
		END FUNCTION funk
	END INTERFACE
	REAL(SP) :: a,b,del
	INTEGER(I4B) :: it
	REAL(SP), DIMENSION(2*3**(n-2)) :: x
	b=exp(-aa)
	a=0.0
	if (n == 1) then
		s=(b-a)*sum(func( (/0.5_sp*(a+b)/) ))
	else
		it=3**(n-2)
		del=(b-a)/(3.0_sp*it)
		x(1:2*it-1:2)=arth(a+0.5_sp*del,3.0_sp*del,it)
		x(2:2*it:2)=x(1:2*it-1:2)+2.0_sp*del
		s=s/3.0_sp+del*sum(func(x))
	end if
	CONTAINS
!BL
		FUNCTION func(x)
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: func
		func=funk(-log(x))/x
		END FUNCTION func
	END SUBROUTINE midexp
