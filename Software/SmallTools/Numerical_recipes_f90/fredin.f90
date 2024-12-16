	FUNCTION fredin(x,a,b,t,f,w,g,ak)
	USE nrtype; USE nrutil, ONLY : assert_eq
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: a,b
	REAL(SP), DIMENSION(:), INTENT(IN) :: x,t,f,w
	REAL(SP), DIMENSION(size(x)) :: fredin
	INTERFACE
		FUNCTION g(t)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), DIMENSION(:), INTENT(IN) :: t
		REAL(SP), DIMENSION(size(t)) :: g
		END FUNCTION g
!BL
		FUNCTION ak(t,s)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), DIMENSION(:), INTENT(IN) :: t,s
		REAL(SP), DIMENSION(size(t),size(s)) :: ak
		END FUNCTION ak
	END INTERFACE
	INTEGER(I4B) :: n
	n=assert_eq(size(f),size(t),size(w),'fredin')
	fredin=g(x)+matmul(ak(x,t),w*f)
	END FUNCTION fredin
