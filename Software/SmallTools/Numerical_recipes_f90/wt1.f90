	SUBROUTINE wt1(a,isign,wtstep)
	USE nrtype; USE nrutil, ONLY : assert
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: a
	INTEGER(I4B), INTENT(IN) :: isign
	INTERFACE
		SUBROUTINE wtstep(a,isign)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: a
		INTEGER(I4B), INTENT(IN) :: isign
		END SUBROUTINE wtstep
	END INTERFACE
	INTEGER(I4B) :: n,nn
	n=size(a)
	call assert(iand(n,n-1)==0, 'n must be a power of 2 in wt1')
	if (n < 4) RETURN
	if (isign >= 0) then
		nn=n
		do
			if (nn < 4) exit
			call wtstep(a(1:nn),isign)
			nn=nn/2
		end do
	else
		nn=4
		do
			if (nn > n) exit
			call wtstep(a(1:nn),isign)
			nn=nn*2
		end do
	end if
	END SUBROUTINE wt1
