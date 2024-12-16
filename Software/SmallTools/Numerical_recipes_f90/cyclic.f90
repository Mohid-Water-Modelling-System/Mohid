	SUBROUTINE cyclic(a,b,c,alpha,beta,r,x)
	USE nrtype; USE nrutil, ONLY : assert,assert_eq
	USE nr, ONLY : tridag
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN):: a,b,c,r
	REAL(SP), INTENT(IN) :: alpha,beta
	REAL(SP), DIMENSION(:), INTENT(OUT):: x
	INTEGER(I4B) :: n
	REAL(SP) :: fact,gamma
	REAL(SP), DIMENSION(size(x)) :: bb,u,z
	n=assert_eq((/size(a),size(b),size(c),size(r),size(x)/),'cyclic')
	call assert(n > 2, 'cyclic arg')
	gamma=-b(1)
	bb(1)=b(1)-gamma
	bb(n)=b(n)-alpha*beta/gamma
	bb(2:n-1)=b(2:n-1)
	call tridag(a(2:n),bb,c(1:n-1),r,x)
	u(1)=gamma
	u(n)=alpha
	u(2:n-1)=0.0
	call tridag(a(2:n),bb,c(1:n-1),u,z)
	fact=(x(1)+beta*x(n)/gamma)/(1.0_sp+z(1)+beta*z(n)/gamma)
	x=x-fact*z
	END SUBROUTINE cyclic
