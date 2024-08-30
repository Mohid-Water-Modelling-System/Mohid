	FUNCTION lop(u)
	USE nrtype; USE nrutil, ONLY : assert_eq
	IMPLICIT NONE
	REAL(DP), DIMENSION(:,:), INTENT(IN) :: u
	REAL(DP), DIMENSION(size(u,1),size(u,1)) :: lop
	INTEGER(I4B) :: n
	REAL(DP) :: h,h2i
	n=assert_eq(size(u,1),size(u,2),'lop')
	h=1.0_dp/(n-1)
	h2i=1.0_dp/(h*h)
	lop(2:n-1,2:n-1)=h2i*(u(3:n,2:n-1)+u(1:n-2,2:n-1)+u(2:n-1,3:n)+&
		u(2:n-1,1:n-2)-4.0_dp*u(2:n-1,2:n-1))+u(2:n-1,2:n-1)**2
	lop(1:n,1)=0.0
	lop(1:n,n)=0.0
	lop(1,1:n)=0.0
	lop(n,1:n)=0.0
	END FUNCTION lop
