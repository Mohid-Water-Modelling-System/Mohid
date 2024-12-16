	FUNCTION bessy_s(n,x)
	USE nrtype; USE nrutil, ONLY : assert
	USE nr, ONLY : bessy0,bessy1
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: n
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: bessy_s
	INTEGER(I4B) :: j
	REAL(SP) :: by,bym,byp,tox
	call assert(n >= 2, x > 0.0, 'bessy_s args')
	tox=2.0_sp/x
	by=bessy1(x)
	bym=bessy0(x)
	do j=1,n-1
		byp=j*tox*by-bym
		bym=by
		by=byp
	end do
	bessy_s=by
	END FUNCTION bessy_s


	FUNCTION bessy_v(n,x)
	USE nrtype; USE nrutil, ONLY : assert
	USE nr, ONLY : bessy0,bessy1
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: n
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(size(x)) :: bessy_v
	INTEGER(I4B) :: j
	REAL(SP), DIMENSION(size(x)) :: by,bym,byp,tox
	call assert(n >= 2, all(x > 0.0), 'bessy_v args')
	tox=2.0_sp/x
	by=bessy1(x)
	bym=bessy0(x)
	do j=1,n-1
		byp=j*tox*by-bym
		bym=by
		by=byp
	end do
	bessy_v=by
	END FUNCTION bessy_v
