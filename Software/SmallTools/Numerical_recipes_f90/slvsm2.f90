	SUBROUTINE slvsm2(u,rhs)
	USE nrtype
	IMPLICIT NONE
	REAL(DP), DIMENSION(3,3), INTENT(OUT) :: u
	REAL(DP), DIMENSION(3,3), INTENT(IN) :: rhs
	REAL(DP) :: disc,fact,h
	u=0.0
	h=0.5_dp
	fact=2.0_dp/h**2
	disc=sqrt(fact**2+rhs(2,2))
	u(2,2)=-rhs(2,2)/(fact+disc)
	END SUBROUTINE slvsm2
