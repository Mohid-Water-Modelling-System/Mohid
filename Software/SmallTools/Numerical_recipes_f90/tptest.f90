	SUBROUTINE tptest(data1,data2,t,prob)
	USE nrtype; USE nrutil, ONLY : assert_eq
	USE nr, ONLY : avevar,betai
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: data1,data2
	REAL(SP), INTENT(OUT) :: t,prob
	INTEGER(I4B) :: n
	REAL(SP) :: ave1,ave2,cov,df,sd,var1,var2
	n=assert_eq(size(data1),size(data2),'tptest')
	call avevar(data1,ave1,var1)
	call avevar(data2,ave2,var2)
	cov=dot_product(data1(:)-ave1,data2(:)-ave2)
	df=n-1
	cov=cov/df
	sd=sqrt((var1+var2-2.0_sp*cov)/n)
	t=(ave1-ave2)/sd
	prob=betai(0.5_sp*df,0.5_sp,df/(df+t**2))
	END SUBROUTINE tptest
