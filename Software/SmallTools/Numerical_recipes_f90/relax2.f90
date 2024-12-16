	SUBROUTINE relax2(u,rhs)
	USE nrtype; USE nrutil, ONLY : assert_eq
	IMPLICIT NONE
	REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: u
	REAL(DP), DIMENSION(:,:), INTENT(IN) :: rhs
	INTEGER(I4B) :: n
	REAL(DP) :: foh2,h,h2i
	REAL(DP) :: res(size(u,1),size(u,1))
	n=assert_eq(size(u,1),size(u,2),size(rhs,1),size(rhs,2),'relax2')
	h=1.0_dp/(n-1)
	h2i=1.0_dp/(h*h)
	foh2=-4.0_dp*h2i
	res(2:n-1:2,2:n-1:2)=h2i*(u(3:n:2,2:n-1:2)+u(1:n-2:2,2:n-1:2)+&
		u(2:n-1:2,3:n:2)+u(2:n-1:2,1:n-2:2)-4.0_dp*u(2:n-1:2,2:n-1:2))&
		+u(2:n-1:2,2:n-1:2)**2-rhs(2:n-1:2,2:n-1:2)
	u(2:n-1:2,2:n-1:2)=u(2:n-1:2,2:n-1:2)-res(2:n-1:2,2:n-1:2)/&
		(foh2+2.0_dp*u(2:n-1:2,2:n-1:2))
	res(3:n-2:2,3:n-2:2)=h2i*(u(4:n-1:2,3:n-2:2)+u(2:n-3:2,3:n-2:2)+&
		u(3:n-2:2,4:n-1:2)+u(3:n-2:2,2:n-3:2)-4.0_dp*u(3:n-2:2,3:n-2:2))&
		+u(3:n-2:2,3:n-2:2)**2-rhs(3:n-2:2,3:n-2:2)
	u(3:n-2:2,3:n-2:2)=u(3:n-2:2,3:n-2:2)-res(3:n-2:2,3:n-2:2)/&
		(foh2+2.0_dp*u(3:n-2:2,3:n-2:2))
	res(3:n-2:2,2:n-1:2)=h2i*(u(4:n-1:2,2:n-1:2)+u(2:n-3:2,2:n-1:2)+&
		u(3:n-2:2,3:n:2)+u(3:n-2:2,1:n-2:2)-4.0_dp*u(3:n-2:2,2:n-1:2))&
		+u(3:n-2:2,2:n-1:2)**2-rhs(3:n-2:2,2:n-1:2)
	u(3:n-2:2,2:n-1:2)=u(3:n-2:2,2:n-1:2)-res(3:n-2:2,2:n-1:2)/&
		(foh2+2.0_dp*u(3:n-2:2,2:n-1:2))
	res(2:n-1:2,3:n-2:2)=h2i*(u(3:n:2,3:n-2:2)+u(1:n-2:2,3:n-2:2)+&
		u(2:n-1:2,4:n-1:2)+u(2:n-1:2,2:n-3:2)-4.0_dp*u(2:n-1:2,3:n-2:2))&
		+u(2:n-1:2,3:n-2:2)**2-rhs(2:n-1:2,3:n-2:2)
	u(2:n-1:2,3:n-2:2)=u(2:n-1:2,3:n-2:2)-res(2:n-1:2,3:n-2:2)/&
		(foh2+2.0_dp*u(2:n-1:2,3:n-2:2))
	END SUBROUTINE relax2
