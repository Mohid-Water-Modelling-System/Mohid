	SUBROUTINE slvsml(u,rhs)
	USE nrtype
	IMPLICIT NONE
	REAL(DP), DIMENSION(3,3), INTENT(OUT) :: u
	REAL(DP), DIMENSION(3,3), INTENT(IN) :: rhs
	REAL(DP) :: h
	u=0.0
	h=0.5_dp
	u(2,2)=-h*h*rhs(2,2)/4.0_dp
	END SUBROUTINE slvsml
