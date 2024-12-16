	FUNCTION gamdev(ia)
	USE nrtype; USE nrutil, ONLY : assert
	USE nr, ONLY : ran1
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: ia
	REAL(SP) :: gamdev
	REAL(SP) :: am,e,h,s,x,y,v(2),arr(5)
	call assert(ia >= 1, 'gamdev arg')
	if (ia < 6) then
		call ran1(arr(1:ia))
		x=-log(product(arr(1:ia)))
	else
		do
			call ran1(v)
			v(2)=2.0_sp*v(2)-1.0_sp
			if (dot_product(v,v) > 1.0) cycle
			y=v(2)/v(1)
			am=ia-1
			s=sqrt(2.0_sp*am+1.0_sp)
			x=s*y+am
			if (x <= 0.0) cycle
			e=(1.0_sp+y**2)*exp(am*log(x/am)-s*y)
			call ran1(h)
			if (h <= e) exit
		end do
	end if
	gamdev=x
	END FUNCTION gamdev
