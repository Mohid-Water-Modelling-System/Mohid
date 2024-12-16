	SUBROUTINE flmoon(n,nph,jd,frac)
	USE nrtype; USE nrutil, ONLY : nrerror
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: n,nph
	INTEGER(I4B), INTENT(OUT) :: jd
	REAL(SP), INTENT(OUT) :: frac
	REAL(SP), PARAMETER :: RAD=PI/180.0_sp
	INTEGER(I4B) :: i
	REAL(SP) :: am,as,c,t,t2,xtra
	c=n+nph/4.0_sp
	t=c/1236.85_sp
	t2=t**2
	as=359.2242_sp+29.105356_sp*c
	am=306.0253_sp+385.816918_sp*c+0.010730_sp*t2
	jd=2415020+28*n+7*nph
	xtra=0.75933_sp+1.53058868_sp*c+(1.178e-4_sp-1.55e-7_sp*t)*t2
	select case(nph)
		case(0,2)
			xtra=xtra+(0.1734_sp-3.93e-4_sp*t)*sin(RAD*as)-0.4068_sp*sin(RAD*am)
		case(1,3)
			xtra=xtra+(0.1721_sp-4.0e-4_sp*t)*sin(RAD*as)-0.6280_sp*sin(RAD*am)
		case default
			call nrerror('flmoon: nph is unknown')
	end select
	i=int(merge(xtra,xtra-1.0_sp, xtra >= 0.0))
	jd=jd+i
	frac=xtra-i
	END SUBROUTINE flmoon
