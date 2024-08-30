	FUNCTION elle_s(phi,ak)
	USE nrtype
	USE nr, ONLY : rd,rf
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: phi,ak
	REAL(SP) :: elle_s
	REAL(SP) :: cc,q,s
	s=sin(phi)
	cc=cos(phi)**2
	q=(1.0_sp-s*ak)*(1.0_sp+s*ak)
	elle_s=s*(rf(cc,q,1.0_sp)-((s*ak)**2)*rd(cc,q,1.0_sp)/3.0_sp)
	END FUNCTION elle_s


	FUNCTION elle_v(phi,ak)
	USE nrtype; USE nrutil, ONLY : assert_eq
	USE nr, ONLY : rd,rf
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: phi,ak
	REAL(SP), DIMENSION(size(phi)) :: elle_v
	REAL(SP), DIMENSION(size(phi)) :: cc,q,s
	INTEGER(I4B) :: ndum
	ndum=assert_eq(size(phi),size(ak),'elle_v')
	s=sin(phi)
	cc=cos(phi)**2
	q=(1.0_sp-s*ak)*(1.0_sp+s*ak)
	elle_v=s*(rf(cc,q,spread(1.0_sp,1,size(phi)))-((s*ak)**2)*&
		rd(cc,q,spread(1.0_sp,1,size(phi)))/3.0_sp)
	END FUNCTION elle_v
