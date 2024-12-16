	SUBROUTINE tutest(data1,data2,t,prob)
	USE nrtype
	USE nr, ONLY : avevar,betai
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: data1,data2
	REAL(SP), INTENT(OUT) :: t,prob
	INTEGER(I4B) :: n1,n2
	REAL(SP) :: ave1,ave2,df,var1,var2
	n1=size(data1)
	n2=size(data2)
	call avevar(data1,ave1,var1)
	call avevar(data2,ave2,var2)
	t=(ave1-ave2)/sqrt(var1/n1+var2/n2)
	df=(var1/n1+var2/n2)**2/((var1/n1)**2/(n1-1)+(var2/n2)**2/(n2-1))
	prob=betai(0.5_sp*df,0.5_sp,df/(df+t**2))
	END SUBROUTINE tutest
