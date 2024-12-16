	FUNCTION resid(u,rhs)
	USE nrtype; USE nrutil, ONLY : assert_eq
	IMPLICIT NONE
	REAL(DP), DIMENSION(:,:), INTENT(IN) :: u,rhs
	REAL(DP), DIMENSION(size(u,1),size(u,1)) :: resid
	INTEGER(I4B) :: n
	REAL(DP) :: h,h2i
	n=assert_eq((/size(u,1),size(u,2),size(rhs,1),size(rhs,2)/),'resid')
	n=size(u,1)
	h=1.0_dp/(n-1)
	h2i=1.0_dp/(h*h)
	resid(2:n-1,2:n-1)=-h2i*(u(3:n,2:n-1)+u(1:n-2,2:n-1)+u(2:n-1,3:n)+&
		u(2:n-1,1:n-2)-4.0_dp*u(2:n-1,2:n-1))+rhs(2:n-1,2:n-1)
	resid(1:n,1)=0.0
	resid(1:n,n)=0.0
	resid(1,1:n)=0.0
	resid(n,1:n)=0.0
	END FUNCTION resid
