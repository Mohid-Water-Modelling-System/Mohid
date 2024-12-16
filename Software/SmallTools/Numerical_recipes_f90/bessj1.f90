	FUNCTION bessj1_s(x)
	USE nrtype; USE nrutil, ONLY : poly
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: bessj1_s
	REAL(SP) :: ax,xx,z
	REAL(DP) :: y
	REAL(DP), DIMENSION(6) :: r = (/72362614232.0_dp,&
		-7895059235.0_dp,242396853.1_dp,-2972611.439_dp,&
		15704.48260_dp,-30.16036606_dp/)
	REAL(DP), DIMENSION(6) :: s = (/144725228442.0_dp,2300535178.0_dp,&
		18583304.74_dp,99447.43394_dp,376.9991397_dp,1.0_dp/)
	REAL(DP), DIMENSION(5) :: p = (/1.0_dp,0.183105e-2_dp,&
		-0.3516396496e-4_dp,0.2457520174e-5_dp,-0.240337019e-6_dp/)
	REAL(DP), DIMENSION(5) :: q = (/0.04687499995_dp,&
		-0.2002690873e-3_dp,0.8449199096e-5_dp,-0.88228987e-6_dp,&
		0.105787412e-6_dp/)
	if (abs(x) < 8.0) then
		y=x**2
		bessj1_s=x*(poly(y,r)/poly(y,s))
	else
		ax=abs(x)
		z=8.0_sp/ax
		y=z**2
		xx=ax-2.356194491_sp
		bessj1_s=sqrt(0.636619772_sp/ax)*(cos(xx)*&
			poly(y,p)-z*sin(xx)*poly(y,q))*sign(1.0_sp,x)
	end if
	END FUNCTION bessj1_s


	FUNCTION bessj1_v(x)
	USE nrtype; USE nrutil, ONLY : poly
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(size(x)) :: bessj1_v
	REAL(SP), DIMENSION(size(x)) :: ax,xx,z
	REAL(DP), DIMENSION(size(x)) :: y
	LOGICAL(LGT), DIMENSION(size(x)) :: mask
	REAL(DP), DIMENSION(6) :: r = (/72362614232.0_dp,&
		-7895059235.0_dp,242396853.1_dp,-2972611.439_dp,&
		15704.48260_dp,-30.16036606_dp/)
	REAL(DP), DIMENSION(6) :: s = (/144725228442.0_dp,2300535178.0_dp,&
		18583304.74_dp,99447.43394_dp,376.9991397_dp,1.0_dp/)
	REAL(DP), DIMENSION(5) :: p = (/1.0_dp,0.183105e-2_dp,&
		-0.3516396496e-4_dp,0.2457520174e-5_dp,-0.240337019e-6_dp/)
	REAL(DP), DIMENSION(5) :: q = (/0.04687499995_dp,&
		-0.2002690873e-3_dp,0.8449199096e-5_dp,-0.88228987e-6_dp,&
		0.105787412e-6_dp/)
	mask = (abs(x) < 8.0)
	where (mask)
		y=x**2
		bessj1_v=x*(poly(y,r,mask)/poly(y,s,mask))
	elsewhere
		ax=abs(x)
		z=8.0_sp/ax
		y=z**2
		xx=ax-2.356194491_sp
		bessj1_v=sqrt(0.636619772_sp/ax)*(cos(xx)*&
			poly(y,p,.not. mask)-z*sin(xx)*poly(y,q,.not. mask))*&
			sign(1.0_sp,x)
	end where
	END FUNCTION bessj1_v
