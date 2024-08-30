	SUBROUTINE ftest(data1,data2,f,prob)
	USE nrtype
	USE nr, ONLY : avevar,betai
	IMPLICIT NONE
	REAL(SP), INTENT(OUT) :: f,prob
	REAL(SP), DIMENSION(:), INTENT(IN) :: data1,data2
	INTEGER(I4B) :: n1,n2
	REAL(SP) :: ave1,ave2,df1,df2,var1,var2
	n1=size(data1)
	n2=size(data2)
	call avevar(data1,ave1,var1)
	call avevar(data2,ave2,var2)
	if (var1 > var2) then
		f=var1/var2
		df1=n1-1
		df2=n2-1
	else
		f=var2/var1
		df1=n2-1
		df2=n1-1
	end if
	prob=2.0_sp*betai(0.5_sp*df2,0.5_sp*df1,df2/(df2+df1*f))
	if (prob > 1.0) prob=2.0_sp-prob
	END SUBROUTINE ftest
