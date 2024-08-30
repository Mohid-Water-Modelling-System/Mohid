	SUBROUTINE kendl1(data1,data2,tau,z,prob)
	USE nrtype; USE nrutil, ONLY : assert_eq
	USE nr, ONLY : erfcc
	IMPLICIT NONE
	REAL(SP), INTENT(OUT) :: tau,z,prob
	REAL(SP), DIMENSION(:), INTENT(IN) :: data1,data2
	INTEGER(I4B) :: is,j,n,n1,n2
	REAL(SP) :: var
	REAL(SP), DIMENSION(size(data1)) :: a1,a2
	n=assert_eq(size(data1),size(data2),'kendl1')
	n1=0
	n2=0
	is=0
	do j=1,n-1
		a1(j+1:n)=data1(j)-data1(j+1:n)
		a2(j+1:n)=data2(j)-data2(j+1:n)
		n1=n1+count(a1(j+1:n) /= 0.0)
		n2=n2+count(a2(j+1:n) /= 0.0)
		is=is+count((a1(j+1:n) > 0.0 .and. a2(j+1:n) > 0.0) &
			.or. (a1(j+1:n) < 0.0 .and. a2(j+1:n) < 0.0)) - &
			count((a1(j+1:n) > 0.0 .and. a2(j+1:n) < 0.0) &
			.or. (a1(j+1:n) < 0.0 .and. a2(j+1:n) > 0.0))
	end do
	tau=real(is,sp)/sqrt(real(n1,sp)*real(n2,sp))
	var=(4.0_sp*n+10.0_sp)/(9.0_sp*n*(n-1.0_sp))
	z=tau/sqrt(var)
	prob=erfcc(abs(z)/SQRT2)
	END SUBROUTINE kendl1
