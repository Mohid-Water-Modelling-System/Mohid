	FUNCTION ellpi_s(phi,en,ak)
	USE nrtype
	USE nr, ONLY : rf,rj
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: phi,en,ak
	REAL(SP) :: ellpi_s
	REAL(SP) :: cc,enss,q,s
	s=sin(phi)
	enss=en*s*s
	cc=cos(phi)**2
	q=(1.0_sp-s*ak)*(1.0_sp+s*ak)
	ellpi_s=s*(rf(cc,q,1.0_sp)-enss*rj(cc,q,1.0_sp,1.0_sp+enss)/3.0_sp)
	END FUNCTION ellpi_s


	FUNCTION ellpi_v(phi,en,ak)
	USE nrtype
	USE nr, ONLY : rf,rj
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: phi,en,ak
	REAL(SP), DIMENSION(size(phi)) :: ellpi_v
	REAL(SP), DIMENSION(size(phi)) :: cc,enss,q,s
	s=sin(phi)
	enss=en*s*s
	cc=cos(phi)**2
	q=(1.0_sp-s*ak)*(1.0_sp+s*ak)
	ellpi_v=s*(rf(cc,q,spread(1.0_sp,1,size(phi)))-enss*&
		rj(cc,q,spread(1.0_sp,1,size(phi)),1.0_sp+enss)/3.0_sp)
	END FUNCTION ellpi_v
