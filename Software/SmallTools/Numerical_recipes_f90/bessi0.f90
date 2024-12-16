	FUNCTION bessi0_s(x)
	USE nrtype; USE nrutil, ONLY : poly
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: bessi0_s
	REAL(SP) :: ax
	REAL(DP), DIMENSION(7) :: p = (/1.0_dp,3.5156229_dp,&
		3.0899424_dp,1.2067492_dp,0.2659732_dp,0.360768e-1_dp,&
		0.45813e-2_dp/)
	REAL(DP), DIMENSION(9) :: q = (/0.39894228_dp,0.1328592e-1_dp,&
		0.225319e-2_dp,-0.157565e-2_dp,0.916281e-2_dp,&
		-0.2057706e-1_dp,0.2635537e-1_dp,-0.1647633e-1_dp,&
		0.392377e-2_dp/)
	ax=abs(x)
	if (ax < 3.75) then
		bessi0_s=poly(real((x/3.75_sp)**2,dp),p)
	else
		bessi0_s=(exp(ax)/sqrt(ax))*poly(real(3.75_sp/ax,dp),q)
	end if
	END FUNCTION bessi0_s


	FUNCTION bessi0_v(x)
	USE nrtype; USE nrutil, ONLY : poly
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(size(x)) :: bessi0_v
	REAL(SP), DIMENSION(size(x)) :: ax
	REAL(DP), DIMENSION(size(x)) :: y
	LOGICAL(LGT), DIMENSION(size(x)) :: mask
	REAL(DP), DIMENSION(7) :: p = (/1.0_dp,3.5156229_dp,&
		3.0899424_dp,1.2067492_dp,0.2659732_dp,0.360768e-1_dp,&
		0.45813e-2_dp/)
	REAL(DP), DIMENSION(9) :: q = (/0.39894228_dp,0.1328592e-1_dp,&
		0.225319e-2_dp,-0.157565e-2_dp,0.916281e-2_dp,&
		-0.2057706e-1_dp,0.2635537e-1_dp,-0.1647633e-1_dp,&
		0.392377e-2_dp/)
	ax=abs(x)
	mask = (ax < 3.75)
	where (mask)
		bessi0_v=poly(real((x/3.75_sp)**2,dp),p,mask)
	elsewhere
		y=3.75_sp/ax
		bessi0_v=(exp(ax)/sqrt(ax))*poly(real(y,dp),q,.not. mask)
	end where
	END FUNCTION bessi0_v
