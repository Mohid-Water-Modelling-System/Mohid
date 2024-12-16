	FUNCTION bessk1_s(x)
	USE nrtype; USE nrutil, ONLY : assert,poly
	USE nr, ONLY : bessi1
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: bessk1_s
	REAL(DP) :: y
	REAL(DP), DIMENSION(7) :: p = (/1.0_dp,0.15443144_dp,&
		-0.67278579_dp,-0.18156897_dp,-0.1919402e-1_dp,&
		-0.110404e-2_dp,-0.4686e-4_dp/)
	REAL(DP), DIMENSION(7) :: q = (/1.25331414_dp,0.23498619_dp,&
		-0.3655620e-1_dp,0.1504268e-1_dp,-0.780353e-2_dp,&
		0.325614e-2_dp,-0.68245e-3_dp/)
	call assert(x > 0.0, 'bessk1_s arg')
	if (x <= 2.0) then
		y=x*x/4.0_sp
		bessk1_s=(log(x/2.0_sp)*bessi1(x))+(1.0_sp/x)*poly(y,p)
	else
		y=2.0_sp/x
		bessk1_s=(exp(-x)/sqrt(x))*poly(y,q)
	end if
	END FUNCTION bessk1_s


	FUNCTION bessk1_v(x)
	USE nrtype; USE nrutil, ONLY : assert,poly
	USE nr, ONLY : bessi1
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(size(x)) :: bessk1_v
	REAL(DP), DIMENSION(size(x)) :: y
	LOGICAL(LGT), DIMENSION(size(x)) :: mask
	REAL(DP), DIMENSION(7) :: p = (/1.0_dp,0.15443144_dp,&
		-0.67278579_dp,-0.18156897_dp,-0.1919402e-1_dp,&
		-0.110404e-2_dp,-0.4686e-4_dp/)
	REAL(DP), DIMENSION(7) :: q = (/1.25331414_dp,0.23498619_dp,&
		-0.3655620e-1_dp,0.1504268e-1_dp,-0.780353e-2_dp,&
		0.325614e-2_dp,-0.68245e-3_dp/)
	call assert(all(x > 0.0), 'bessk1_v arg')
	mask = (x <= 2.0)
	where (mask)
		y=x*x/4.0_sp
		bessk1_v=(log(x/2.0_sp)*bessi1(x))+(1.0_sp/x)*poly(y,p,mask)
	elsewhere
		y=2.0_sp/x
		bessk1_v=(exp(-x)/sqrt(x))*poly(y,q,.not. mask)
	end where
	END FUNCTION bessk1_v
