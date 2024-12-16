	RECURSIVE FUNCTION recur1(a,b) RESULT(u)
	USE nrtype; USE nrutil, ONLY : assert_eq
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: a,b
	REAL(SP), DIMENSION(size(a)) :: u
	INTEGER(I4B), PARAMETER :: NPAR_RECUR1=8
	INTEGER(I4B) :: n,j
	n=assert_eq(size(a),size(b)+1,'recur1')
	u(1)=a(1)
	if (n < NPAR_RECUR1) then
		do j=2,n
			u(j)=a(j)+b(j-1)*u(j-1)
		end do
	else
		u(2:n:2)=recur1(a(2:n:2)+a(1:n-1:2)*b(1:n-1:2), &
				b(3:n-1:2)*b(2:n-2:2))
		u(3:n:2)=a(3:n:2)+b(2:n-1:2)*u(2:n-1:2)
	end if
	END FUNCTION recur1
