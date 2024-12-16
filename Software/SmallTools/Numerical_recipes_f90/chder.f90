	FUNCTION chder(a,b,c)
	USE nrtype; USE nrutil, ONLY : arth,cumsum
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: a,b
	REAL(SP), DIMENSION(:), INTENT(IN) :: c
	REAL(SP), DIMENSION(size(c)) :: chder
	INTEGER(I4B) :: n
	REAL(SP) :: con
	REAL(SP), DIMENSION(size(c)) :: temp
	n=size(c)
	temp(1)=0.0
	temp(2:n)=2.0_sp*arth(n-1,-1,n-1)*c(n:2:-1)
	chder(n:1:-2)=cumsum(temp(1:n:2))
	chder(n-1:1:-2)=cumsum(temp(2:n:2))
	con=2.0_sp/(b-a)
	chder=chder*con
	END FUNCTION chder
