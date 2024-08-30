	FUNCTION evlmem(fdt,d,xms)
	USE nrtype; USE nrutil, ONLY : poly
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: fdt,xms
	REAL(SP), DIMENSION(:), INTENT(IN) :: d
	REAL(SP) :: evlmem
	COMPLEX(SPC) :: z,zz
	REAL(DP) :: theta
	theta=TWOPI_D*fdt
	z=cmplx(cos(theta),sin(theta),kind=spc)
	zz=1.0_sp-z*poly(z,d)
	evlmem=xms/abs(zz)**2
	END FUNCTION evlmem
