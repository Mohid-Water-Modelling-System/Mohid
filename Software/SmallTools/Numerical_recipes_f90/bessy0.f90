	FUNCTION bessy0_s(x)
	USE nrtype; USE nrutil, ONLY : assert,poly
	USE nr, ONLY : bessj0
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: bessy0_s
	REAL(SP) :: xx,z
	REAL(DP) :: y
	REAL(DP), DIMENSION(5) :: p = (/1.0_dp,-0.1098628627e-2_dp,&
		0.2734510407e-4_dp,-0.2073370639e-5_dp,0.2093887211e-6_dp/)
	REAL(DP), DIMENSION(5) :: q = (/-0.1562499995e-1_dp,&
		0.1430488765e-3_dp,-0.6911147651e-5_dp,0.7621095161e-6_dp,&
		-0.934945152e-7_dp/)
	REAL(DP), DIMENSION(6) :: r = (/-2957821389.0_dp,7062834065.0_dp,&
		-512359803.6_dp,10879881.29_dp,-86327.92757_dp,&
		228.4622733_dp/)
	REAL(DP), DIMENSION(6) :: s = (/40076544269.0_dp,745249964.8_dp,&
		7189466.438_dp,47447.26470_dp,226.1030244_dp,1.0_dp/)
	call assert(x > 0.0, 'bessy0_s arg')
	if (abs(x) < 8.0) then
		y=x**2
		bessy0_s=(poly(y,r)/poly(y,s))+&
			0.636619772_sp*bessj0(x)*log(x)
	else
		z=8.0_sp/x
		y=z**2
		xx=x-0.785398164_sp
		bessy0_s=sqrt(0.636619772_sp/x)*(sin(xx)*&
			poly(y,p)+z*cos(xx)*poly(y,q))
	end if
	END FUNCTION bessy0_s


	FUNCTION bessy0_v(x)
	USE nrtype; USE nrutil, ONLY : assert,poly
	USE nr, ONLY : bessj0
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(size(x)) :: bessy0_v
	REAL(SP), DIMENSION(size(x)) :: xx,z
	REAL(DP), DIMENSION(size(x)) :: y
	LOGICAL(LGT), DIMENSION(size(x)) :: mask
	REAL(DP), DIMENSION(5) :: p = (/1.0_dp,-0.1098628627e-2_dp,&
		0.2734510407e-4_dp,-0.2073370639e-5_dp,0.2093887211e-6_dp/)
	REAL(DP), DIMENSION(5) :: q = (/-0.1562499995e-1_dp,&
		0.1430488765e-3_dp,-0.6911147651e-5_dp,0.7621095161e-6_dp,&
		-0.934945152e-7_dp/)
	REAL(DP), DIMENSION(6) :: r = (/-2957821389.0_dp,7062834065.0_dp,&
		-512359803.6_dp,10879881.29_dp,-86327.92757_dp,&
		228.4622733_dp/)
	REAL(DP), DIMENSION(6) :: s = (/40076544269.0_dp,745249964.8_dp,&
		7189466.438_dp,47447.26470_dp,226.1030244_dp,1.0_dp/)
	call assert(all(x > 0.0), 'bessy0_v arg')
	mask = (abs(x) < 8.0)
	where (mask)
		y=x**2
		bessy0_v=(poly(y,r,mask)/poly(y,s,mask))+&
			0.636619772_sp*bessj0(x)*log(x)
	elsewhere
		z=8.0_sp/x
		y=z**2
		xx=x-0.785398164_sp
		bessy0_v=sqrt(0.636619772_sp/x)*(sin(xx)*&
			poly(y,p,.not. mask)+z*cos(xx)*poly(y,q,.not. mask))
	end where
	END FUNCTION bessy0_v
