MODULE kermom_info
	USE nrtype
	REAL(DP) :: kermom_x
END MODULE kermom_info

	FUNCTION kermom(y,m)
	USE nrtype
	USE kermom_info
	IMPLICIT NONE
	REAL(DP),  INTENT(IN) :: y
	INTEGER(I4B), INTENT(IN) :: m
	REAL(DP), DIMENSION(m) :: kermom
	REAL(DP) :: x,d,df,clog,x2,x3,x4
	x=kermom_x
	if (y >= x) then
		d=y-x
		df=2.0_dp*sqrt(d)*d
		kermom(1:4) = (/ df/3.0_dp, df*(x/3.0_dp+d/5.0_dp),&
			df*((x/3.0_dp + 0.4_dp*d)*x + d**2/7.0_dp),&
			df*(((x/3.0_dp + 0.6_dp*d)*x + 3.0_dp*d**2/7.0_dp)*x&
				+ d**3/9.0_dp) /)
	else
		x2=x**2
		x3=x2*x
		x4=x2*x2
		d=x-y
		clog=log(d)
		kermom(1:4) = (/ d*(clog-1.0_dp),&
			-0.25_dp*(3.0_dp*x+y-2.0_dp*clog*(x+y))*d,&
			(-11.0_dp*x3+y*(6.0_dp*x2+y*(3.0_dp*x+2.0_dp*y))&
				+6.0_dp*clog*(x3-y**3))/18.0_dp,&
			(-25.0_dp*x4+y*(12.0_dp*x3+y*(6.0_dp*x2+y*&
				(4.0_dp*x+3.0_dp*y)))+12.0_dp*clog*(x4-y**4))/48.0_dp /)
	end if
	END FUNCTION kermom
