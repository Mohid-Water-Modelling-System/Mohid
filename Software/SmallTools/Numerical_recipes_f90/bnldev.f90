	FUNCTION bnldev(pp,n)
	USE nrtype
	USE nr, ONLY : gammln,ran1
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: pp
	INTEGER(I4B), INTENT(IN) :: n
	REAL(SP) :: bnldev
	INTEGER(I4B) :: j
	INTEGER(I4B), SAVE :: nold=-1
	REAL(SP) :: am,em,g,h,p,sq,t,y,arr(24)
	REAL(SP), SAVE :: pc,plog,pclog,en,oldg,pold=-1.0
	p=merge(pp,1.0_sp-pp, pp <= 0.5_sp )
	am=n*p
	if (n < 25) then
		call ran1(arr(1:n))
		bnldev=count(arr(1:n)<p)
	else if (am < 1.0) then
		g=exp(-am)
		t=1.0
		do j=0,n
			call ran1(h)
			t=t*h
			if (t < g) exit
		end do
		bnldev=merge(j,n, j <= n)
	else
		if (n /= nold) then
			en=n
			oldg=gammln(en+1.0_sp)
			nold=n
		end if
		if (p /= pold) then
			pc=1.0_sp-p
			plog=log(p)
			pclog=log(pc)
			pold=p
		end if
		sq=sqrt(2.0_sp*am*pc)
		do
			call ran1(h)
			y=tan(PI*h)
			em=sq*y+am
			if (em < 0.0 .or. em >= en+1.0_sp) cycle
			em=int(em)
			t=1.2_sp*sq*(1.0_sp+y**2)*exp(oldg-gammln(em+1.0_sp)-&
				gammln(en-em+1.0_sp)+em*plog+(en-em)*pclog)
			call ran1(h)
			if (h <= t) exit
		end do
		bnldev=em
	end if
	if (p /= pp) bnldev=n-bnldev
	END FUNCTION bnldev
