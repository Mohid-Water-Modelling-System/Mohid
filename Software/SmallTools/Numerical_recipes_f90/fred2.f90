	SUBROUTINE fred2(a,b,t,f,w,g,ak)
	USE nrtype; USE nrutil, ONLY : assert_eq,unit_matrix
	USE nr, ONLY : gauleg,lubksb,ludcmp
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: a,b
	REAL(SP), DIMENSION(:), INTENT(OUT) :: t,f,w
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
	INTEGER(I4B), DIMENSION(size(f)) :: indx
	REAL(SP) :: d
	REAL(SP), DIMENSION(size(f),size(f)) :: omk
	n=assert_eq(size(f),size(t),size(w),'fred2')
	call gauleg(a,b,t,w)
	call unit_matrix(omk)
	omk=omk-ak(t,t)*spread(w,dim=1,ncopies=n)
	f=g(t)
	call ludcmp(omk,indx,d)
	call lubksb(omk,indx,f)
	END SUBROUTINE fred2
