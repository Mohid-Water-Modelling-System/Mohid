	SUBROUTINE mp2dfr(a,s,n,m)
	USE nrtype
	USE mpops, ONLY : mplsh,mpsmu
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: n
	INTEGER(I4B), INTENT(OUT) :: m
	CHARACTER(1), DIMENSION(:), INTENT(INOUT) :: a
	CHARACTER(1), DIMENSION(:), INTENT(OUT) :: s
	INTEGER(I4B), PARAMETER :: IAZ=48
	INTEGER(I4B) :: j
	m=int(2.408_sp*n)
	do j=1,m
		call mpsmu(a,a,n,10)
		s(j)=char(ichar(a(1))+IAZ)
		call mplsh(a,n)
	end do
	END SUBROUTINE mp2dfr
