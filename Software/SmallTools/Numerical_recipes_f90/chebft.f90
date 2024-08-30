	FUNCTION chebft(a,b,n,func)
	USE nrtype; USE nrutil, ONLY : arth,outerprod
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: a,b
	INTEGER(I4B), INTENT(IN) :: n
	REAL(SP), DIMENSION(n) :: chebft
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: func
		END FUNCTION func
	END INTERFACE
	REAL(DP) :: bma,bpa
	REAL(DP), DIMENSION(n) :: theta
	bma=0.5_dp*(b-a)
	bpa=0.5_dp*(b+a)
	theta(:)=PI_D*arth(0.5_dp,1.0_dp,n)/n
	chebft(:)=matmul(cos(outerprod(arth(0.0_dp,1.0_dp,n),theta)), &
		func(real(cos(theta)*bma+bpa,sp)))*2.0_dp/n
	END FUNCTION chebft
