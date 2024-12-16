	SUBROUTINE splie2(x1a,x2a,ya,y2a)
	USE nrtype; USE nrutil, ONLY : assert_eq
	USE nr, ONLY : spline
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x1a,x2a
	REAL(SP), DIMENSION(:,:), INTENT(IN) :: ya
	REAL(SP), DIMENSION(:,:), INTENT(OUT) :: y2a
	INTEGER(I4B) :: j,m,ndum
	m=assert_eq(size(x1a),size(ya,1),size(y2a,1),'splie2: m')
	ndum=assert_eq(size(x2a),size(ya,2),size(y2a,2),'splie2: ndum')
	do j=1,m
		call spline(x2a,ya(j,:),1.0e30_sp,1.0e30_sp,y2a(j,:))
	end do
	END SUBROUTINE splie2
