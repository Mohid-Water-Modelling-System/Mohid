	FUNCTION bessy1_s(x)
	USE nrtype; USE nrutil, ONLY : assert,poly
	USE nr, ONLY : bessj1
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: bessy1_s
	REAL(SP) :: xx,z
	REAL(DP) :: y
	REAL(DP), DIMENSION(5) :: p = (/1.0_dp,0.183105e-2_dp,&
		-0.3516396496e-4_dp,0.2457520174e-5_dp,-0.240337019e-6_dp/)
	REAL(DP), DIMENSION(5) :: q = (/0.04687499995_dp,&
		-0.2002690873e-3_dp,0.8449199096e-5_dp,-0.88228987e-6_dp,&
		0.105787412e-6_dp/)
	REAL(DP), DIMENSION(6) :: r = (/-0.4900604943e13_dp,&
		0.1275274390e13_dp,-0.5153438139e11_dp,0.7349264551e9_dp,&
		-0.4237922726e7_dp,0.8511937935e4_dp/)
	REAL(DP), DIMENSION(7) :: s = (/0.2499580570e14_dp,&
		0.4244419664e12_dp,0.3733650367e10_dp,0.2245904002e8_dp,&
		0.1020426050e6_dp,0.3549632885e3_dp,1.0_dp/)
	call assert(x > 0.0, 'bessy1_s arg')
	if (abs(x) < 8.0) then
		y=x**2
		bessy1_s=x*(poly(y,r)/poly(y,s))+&
			0.636619772_sp*(bessj1(x)*log(x)-1.0_sp/x)
	else
		z=8.0_sp/x
		y=z**2
		xx=x-2.356194491_sp
		bessy1_s=sqrt(0.636619772_sp/x)*(sin(xx)*&
			poly(y,p)+z*cos(xx)*poly(y,q))
	end if
	END FUNCTION bessy1_s


	FUNCTION bessy1_v(x)
	USE nrtype; USE nrutil, ONLY : assert,poly
	USE nr, ONLY : bessj1
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(size(x)) :: bessy1_v
	REAL(SP), DIMENSION(size(x)) :: xx,z
	REAL(DP), DIMENSION(size(x)) :: y
	LOGICAL(LGT), DIMENSION(size(x)) :: mask
	REAL(DP), DIMENSION(5) :: p = (/1.0_dp,0.183105e-2_dp,&
		-0.3516396496e-4_dp,0.2457520174e-5_dp,-0.240337019e-6_dp/)
	REAL(DP), DIMENSION(5) :: q = (/0.04687499995_dp,&
		-0.2002690873e-3_dp,0.8449199096e-5_dp,-0.88228987e-6_dp,&
		0.105787412e-6_dp/)
	REAL(DP), DIMENSION(6) :: r = (/-0.4900604943e13_dp,&
		0.1275274390e13_dp,-0.5153438139e11_dp,0.7349264551e9_dp,&
		-0.4237922726e7_dp,0.8511937935e4_dp/)
	REAL(DP), DIMENSION(7) :: s = (/0.2499580570e14_dp,&
		0.4244419664e12_dp,0.3733650367e10_dp,0.2245904002e8_dp,&
		0.1020426050e6_dp,0.3549632885e3_dp,1.0_dp/)
	call assert(all(x > 0.0), 'bessy1_v arg')
	mask = (abs(x) < 8.0)
	where (mask)
		y=x**2
		bessy1_v=x*(poly(y,r,mask)/poly(y,s,mask))+&
			0.636619772_sp*(bessj1(x)*log(x)-1.0_sp/x)
	elsewhere
		z=8.0_sp/x
		y=z**2
		xx=x-2.356194491_sp
		bessy1_v=sqrt(0.636619772_sp/x)*(sin(xx)*&
			poly(y,p,.not. mask)+z*cos(xx)*poly(y,q,.not. mask))
	end where
	END FUNCTION bessy1_v
