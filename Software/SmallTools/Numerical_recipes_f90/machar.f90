	SUBROUTINE machar(ibeta,it,irnd,ngrd,machep,negep,iexp,minexp,&
		maxexp,eps,epsneg,xmin,xmax)
	USE nrtype
	IMPLICIT NONE
	INTEGER(I4B), INTENT(OUT) :: ibeta,iexp,irnd,it,machep,maxexp,minexp,negep,ngrd
	REAL(SP), INTENT(OUT) :: eps,epsneg,xmax,xmin
	REAL(SP), PARAMETER :: RX=1.0
	REAL(SP) :: a,beta,betah,one,temp,tempa,two,zero
	ibeta=radix(RX)
	it=digits(RX)
	machep=exponent(nearest(RX,RX)-RX)-1
	negep=exponent(nearest(RX,-RX)-RX)-1
	minexp=minexponent(RX)-1
	maxexp=maxexponent(RX)
	iexp=nint(log(real(maxexp-minexp+2,sp))/log(2.0_sp))
	eps=real(ibeta,sp)**machep
	epsneg=real(ibeta,sp)**negep
	xmax=huge(RX)
	xmin=tiny(RX)
	one=RX
	two=one+one
	zero=one-one
	beta=real(ibeta,sp)
	a=beta**(-negep)
	irnd=0
	betah=beta/two
	temp=a+betah
	if (temp-a /= zero) irnd=1
	tempa=a+beta
	temp=tempa+betah
	if ((irnd == 0) .and. (temp-tempa /= zero)) irnd=2
	ngrd=0
	temp=one+eps
	if ((irnd == 0) .and. (temp*one-one /= zero)) ngrd=1
	temp=xmin/two
	if (temp /= zero) irnd=irnd+3
	END SUBROUTINE machar
