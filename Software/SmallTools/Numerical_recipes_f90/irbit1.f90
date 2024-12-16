	FUNCTION irbit1(iseed)
	USE nrtype
	IMPLICIT NONE
	INTEGER(I4B), INTENT(INOUT) :: iseed
	INTEGER(I4B) :: irbit1
	if (btest(iseed,17) .neqv. btest(iseed,4) .neqv. btest(iseed,1) &
		.neqv. btest(iseed,0)) then
		iseed=ibset(ishft(iseed,1),0)
		irbit1=1
	else
		iseed=ishft(iseed,1)
		irbit1=0
	end if
	END FUNCTION irbit1
