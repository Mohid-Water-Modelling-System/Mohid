	FUNCTION chint(a,b,c)
	USE nrtype; USE nrutil, ONLY : arth
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: a,b
	REAL(SP), DIMENSION(:), INTENT(IN) :: c
	REAL(SP), DIMENSION(size(c)) :: chint
	INTEGER(I4B) :: n
	REAL(SP) :: con
	n=size(c)
	con=0.25_sp*(b-a)
	chint(2:n-1)=con*(c(1:n-2)-c(3:n))/arth(1,1,n-2)
	chint(n)=con*c(n-1)/(n-1)
	chint(1)=2.0_sp*(sum(chint(2:n:2))-sum(chint(3:n:2)))
	END FUNCTION chint
