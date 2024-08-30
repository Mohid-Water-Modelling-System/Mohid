	SUBROUTINE pcshft(a,b,d)
	USE nrtype; USE nrutil, ONLY : geop
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: a,b
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: d
	INTEGER(I4B) :: j,n
	REAL(SP), DIMENSION(size(d)) :: dd
	REAL(SP) :: x
	n=size(d)
	dd=d*geop(1.0_sp,2.0_sp/(b-a),n)
	x=-0.5_sp*(a+b)
	d(1)=dd(n)
	d(2:n)=0.0
	do j=n-1,1,-1
		d(2:n+1-j)=d(2:n+1-j)*x+d(1:n-j)
		d(1)=d(1)*x+dd(j)
	end do
	END SUBROUTINE pcshft
