	FUNCTION igray(n,is)
	USE nrtype
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: n,is
	INTEGER(I4B) :: igray
	INTEGER(I4B) :: idiv,ish
	if (is >= 0) then
		igray=ieor(n,n/2)
	else
		ish=-1
		igray=n
		do
			idiv=ishft(igray,ish)
			igray=ieor(igray,idiv)
			if (idiv <= 1 .or. ish == -16) RETURN
			ish=ish+ish
		end do
	end if
	END FUNCTION igray
