	FUNCTION snrm(sx,itol)
	USE nrtype
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), INTENT(IN) :: sx
	INTEGER(I4B), INTENT(IN) :: itol
	REAL(DP) :: snrm
	if (itol <= 3) then
		snrm=sqrt(dot_product(sx,sx))
	else
		snrm=maxval(abs(sx))
	end if
	END FUNCTION snrm
