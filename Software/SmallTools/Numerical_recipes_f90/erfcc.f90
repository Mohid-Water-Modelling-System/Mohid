	FUNCTION erfcc_s(x)
	USE nrtype; USE nrutil, ONLY : poly
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: erfcc_s
	REAL(SP) :: t,z
	REAL(SP), DIMENSION(10) :: coef = (/-1.26551223_sp,1.00002368_sp,&
		0.37409196_sp,0.09678418_sp,-0.18628806_sp,0.27886807_sp,&
		-1.13520398_sp,1.48851587_sp,-0.82215223_sp,0.17087277_sp/)
	z=abs(x)
	t=1.0_sp/(1.0_sp+0.5_sp*z)
	erfcc_s=t*exp(-z*z+poly(t,coef))
	if (x < 0.0) erfcc_s=2.0_sp-erfcc_s
	END FUNCTION erfcc_s


	FUNCTION erfcc_v(x)
	USE nrtype; USE nrutil, ONLY : poly
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(size(x)) :: erfcc_v,t,z
	REAL(SP), DIMENSION(10) :: coef = (/-1.26551223_sp,1.00002368_sp,&
		0.37409196_sp,0.09678418_sp,-0.18628806_sp,0.27886807_sp,&
		-1.13520398_sp,1.48851587_sp,-0.82215223_sp,0.17087277_sp/)
	z=abs(x)
	t=1.0_sp/(1.0_sp+0.5_sp*z)
	erfcc_v=t*exp(-z*z+poly(t,coef))
	where (x < 0.0) erfcc_v=2.0_sp-erfcc_v
	END FUNCTION erfcc_v
