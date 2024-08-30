	SUBROUTINE banmul(a,m1,m2,x,b)
	USE nrtype; USE nrutil, ONLY : assert_eq,arth
	IMPLICIT NONE
	REAL(SP), DIMENSION(:,:), INTENT(IN) :: a
	INTEGER(I4B), INTENT(IN) :: m1,m2
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(:), INTENT(OUT) :: b
	INTEGER(I4B) :: m,n
	n=assert_eq(size(a,1),size(b),size(x),'banmul: n')
	m=assert_eq(size(a,2),m1+m2+1,'banmul: m')
	b=sum(a*eoshift(spread(x,dim=2,ncopies=m), &
		dim=1,shift=arth(-m1,1,m)),dim=2)
	END SUBROUTINE banmul
