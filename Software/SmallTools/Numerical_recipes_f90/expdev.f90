	SUBROUTINE expdev_s(harvest)
	USE nrtype
	USE nr, ONLY : ran1
	IMPLICIT NONE
	REAL(SP), INTENT(OUT) :: harvest
	REAL(SP) :: dum
	call ran1(dum)
	harvest=-log(dum)
	END SUBROUTINE expdev_s

	SUBROUTINE expdev_v(harvest)
	USE nrtype
	USE nr, ONLY : ran1
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(OUT) :: harvest
	REAL(SP), DIMENSION(size(harvest)) :: dum
	call ran1(dum)
	harvest=-log(dum)
	END SUBROUTINE expdev_v
