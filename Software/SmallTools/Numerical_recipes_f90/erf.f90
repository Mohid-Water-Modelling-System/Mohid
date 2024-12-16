	FUNCTION erf_s(x)
	USE nrtype
	USE nr, ONLY : gammp
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: erf_s
	erf_s=gammp(0.5_sp,x**2)
	if (x < 0.0) erf_s=-erf_s
	END FUNCTION erf_s


	FUNCTION erf_v(x)
	USE nrtype
	USE nr, ONLY : gammp
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(size(x)) :: erf_v
	erf_v=gammp(spread(0.5_sp,1,size(x)),x**2)
	where (x < 0.0) erf_v=-erf_v
	END FUNCTION erf_v
