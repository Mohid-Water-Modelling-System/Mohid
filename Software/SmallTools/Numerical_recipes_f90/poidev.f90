	FUNCTION poidev(xm)
	USE nrtype
	USE nr, ONLY : gammln,ran1
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: xm
	REAL(SP) :: poidev
	REAL(SP) :: em,harvest,t,y
	REAL(SP), SAVE :: alxm,g,oldm=-1.0_sp,sq
	if (xm < 12.0) then
		if (xm /= oldm) then
			oldm=xm
			g=exp(-xm)
		end if
		em=-1
		t=1.0
		do
			em=em+1.0_sp
			call ran1(harvest)
			t=t*harvest
			if (t <= g) exit
		end do
	else
		if (xm /= oldm) then
			oldm=xm
			sq=sqrt(2.0_sp*xm)
			alxm=log(xm)
			g=xm*alxm-gammln(xm+1.0_sp)
		end if
		do
			do
				call ran1(harvest)
				y=tan(PI*harvest)
				em=sq*y+xm
				if (em >= 0.0) exit
			end do
			em=int(em)
			t=0.9_sp*(1.0_sp+y**2)*exp(em*alxm-gammln(em+1.0_sp)-g)
			call ran1(harvest)
			if (harvest <= t) exit
		end do
	end if
	poidev=em
	END FUNCTION poidev
