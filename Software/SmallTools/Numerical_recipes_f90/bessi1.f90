	FUNCTION bessi1_s(x)
	USE nrtype; USE nrutil, ONLY : poly
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: bessi1_s
	REAL(SP) :: ax
	REAL(DP), DIMENSION(7) :: p = (/0.5_dp,0.87890594_dp,&
		0.51498869_dp,0.15084934_dp,0.2658733e-1_dp,&
		0.301532e-2_dp,0.32411e-3_dp/)
	REAL(DP), DIMENSION(9) :: q = (/0.39894228_dp,-0.3988024e-1_dp,&
		-0.362018e-2_dp,0.163801e-2_dp,-0.1031555e-1_dp,&
		0.2282967e-1_dp,-0.2895312e-1_dp,0.1787654e-1_dp,&
		-0.420059e-2_dp/)
	ax=abs(x)
	if (ax < 3.75) then
		bessi1_s=ax*poly(real((x/3.75_sp)**2,dp),p)
	else
		bessi1_s=(exp(ax)/sqrt(ax))*poly(real(3.75_sp/ax,dp),q)
	end if
	if (x < 0.0) bessi1_s=-bessi1_s
	END FUNCTION bessi1_s


	FUNCTION bessi1_v(x)
	USE nrtype; USE nrutil, ONLY : poly
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(size(x)) :: bessi1_v
	REAL(SP), DIMENSION(size(x)) :: ax
	REAL(DP), DIMENSION(size(x)) :: y
	LOGICAL(LGT), DIMENSION(size(x)) :: mask
	REAL(DP), DIMENSION(7) :: p = (/0.5_dp,0.87890594_dp,&
		0.51498869_dp,0.15084934_dp,0.2658733e-1_dp,&
		0.301532e-2_dp,0.32411e-3_dp/)
	REAL(DP), DIMENSION(9) :: q = (/0.39894228_dp,-0.3988024e-1_dp,&
		-0.362018e-2_dp,0.163801e-2_dp,-0.1031555e-1_dp,&
		0.2282967e-1_dp,-0.2895312e-1_dp,0.1787654e-1_dp,&
		-0.420059e-2_dp/)
	ax=abs(x)
	mask = (ax < 3.75)
	where (mask)
		bessi1_v=ax*poly(real((x/3.75_sp)**2,dp),p,mask)
	elsewhere
		y=3.75_sp/ax
		bessi1_v=(exp(ax)/sqrt(ax))*poly(real(y,dp),q,.not. mask)
	end where
	where (x < 0.0) bessi1_v=-bessi1_v
	END FUNCTION bessi1_v
