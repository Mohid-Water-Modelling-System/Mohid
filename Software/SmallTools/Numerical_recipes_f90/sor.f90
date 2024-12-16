	SUBROUTINE sor(a,b,c,d,e,f,u,rjac)
	USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
	IMPLICIT NONE
	REAL(DP), DIMENSION(:,:), INTENT(IN) :: a,b,c,d,e,f
	REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: u
	REAL(DP), INTENT(IN) :: rjac
	INTEGER(I4B), PARAMETER :: MAXITS=1000
	REAL(DP), PARAMETER :: EPS=1.0e-5_dp
	REAL(DP), DIMENSION(size(a,1),size(a,1)) :: resid
	INTEGER(I4B) :: jmax,jm1,jm2,jm3,n
	REAL(DP) :: anorm,anormf,omega
	jmax=assert_eq((/size(a,1),size(a,2),size(b,1),size(b,2), &
		size(c,1),size(c,2),size(d,1),size(d,2),size(e,1), &
		size(e,2),size(f,1),size(f,2),size(u,1),size(u,2)/),'sor')
	jm1=jmax-1
	jm2=jmax-2
	jm3=jmax-3
	anormf=sum(abs(f(2:jm1,2:jm1)))
	omega=1.0
	do n=1,MAXITS
		resid(2:jm1:2,2:jm1:2)=a(2:jm1:2,2:jm1:2)*u(3:jmax:2,2:jm1:2)+&
			b(2:jm1:2,2:jm1:2)*u(1:jm2:2,2:jm1:2)+&
			c(2:jm1:2,2:jm1:2)*u(2:jm1:2,3:jmax:2)+&
			d(2:jm1:2,2:jm1:2)*u(2:jm1:2,1:jm2:2)+&
			e(2:jm1:2,2:jm1:2)*u(2:jm1:2,2:jm1:2)-f(2:jm1:2,2:jm1:2)
		u(2:jm1:2,2:jm1:2)=u(2:jm1:2,2:jm1:2)-omega*&
			resid(2:jm1:2,2:jm1:2)/e(2:jm1:2,2:jm1:2)
		resid(3:jm2:2,3:jm2:2)=a(3:jm2:2,3:jm2:2)*u(4:jm1:2,3:jm2:2)+&
			b(3:jm2:2,3:jm2:2)*u(2:jm3:2,3:jm2:2)+&
			c(3:jm2:2,3:jm2:2)*u(3:jm2:2,4:jm1:2)+&
			d(3:jm2:2,3:jm2:2)*u(3:jm2:2,2:jm3:2)+&
			e(3:jm2:2,3:jm2:2)*u(3:jm2:2,3:jm2:2)-f(3:jm2:2,3:jm2:2)
		u(3:jm2:2,3:jm2:2)=u(3:jm2:2,3:jm2:2)-omega*&
			resid(3:jm2:2,3:jm2:2)/e(3:jm2:2,3:jm2:2)
		omega=merge(1.0_dp/(1.0_dp-0.5_dp*rjac**2), &
			1.0_dp/(1.0_dp-0.25_dp*rjac**2*omega), n == 1)
		resid(3:jm2:2,2:jm1:2)=a(3:jm2:2,2:jm1:2)*u(4:jm1:2,2:jm1:2)+&
			b(3:jm2:2,2:jm1:2)*u(2:jm3:2,2:jm1:2)+&
			c(3:jm2:2,2:jm1:2)*u(3:jm2:2,3:jmax:2)+&
			d(3:jm2:2,2:jm1:2)*u(3:jm2:2,1:jm2:2)+&
			e(3:jm2:2,2:jm1:2)*u(3:jm2:2,2:jm1:2)-f(3:jm2:2,2:jm1:2)
		u(3:jm2:2,2:jm1:2)=u(3:jm2:2,2:jm1:2)-omega*&
			resid(3:jm2:2,2:jm1:2)/e(3:jm2:2,2:jm1:2)
		resid(2:jm1:2,3:jm2:2)=a(2:jm1:2,3:jm2:2)*u(3:jmax:2,3:jm2:2)+&
			b(2:jm1:2,3:jm2:2)*u(1:jm2:2,3:jm2:2)+&
			c(2:jm1:2,3:jm2:2)*u(2:jm1:2,4:jm1:2)+&
			d(2:jm1:2,3:jm2:2)*u(2:jm1:2,2:jm3:2)+&
			e(2:jm1:2,3:jm2:2)*u(2:jm1:2,3:jm2:2)-f(2:jm1:2,3:jm2:2)
		u(2:jm1:2,3:jm2:2)=u(2:jm1:2,3:jm2:2)-omega*&
			resid(2:jm1:2,3:jm2:2)/e(2:jm1:2,3:jm2:2)
		omega=1.0_dp/(1.0_dp-0.25_dp*rjac**2*omega)
		anorm=sum(abs(resid(2:jm1,2:jm1)))
		if (anorm < EPS*anormf) exit
	end do
	if (n > MAXITS) call nrerror('MAXITS exceeded in sor')
	END SUBROUTINE sor
