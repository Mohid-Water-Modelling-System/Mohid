	FUNCTION irbit2(iseed)
	USE nrtype
	IMPLICIT NONE
	INTEGER(I4B), INTENT(INOUT) :: iseed
	INTEGER(I4B) :: irbit2
	INTEGER(I4B), PARAMETER :: IB1=1,IB2=2,IB5=16,MASK=IB1+IB2+IB5
	if (btest(iseed,17)) then
		iseed=ibset(ishft(ieor(iseed,MASK),1),0)
		irbit2=1
	else
		iseed=ibclr(ishft(iseed,1),0)
		irbit2=0
	end if
	END FUNCTION irbit2
