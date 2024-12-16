	SUBROUTINE sprstp(sa)
	USE nrtype
	IMPLICIT NONE
	TYPE(sprs2_sp), INTENT(INOUT) :: sa
	INTEGER(I4B), DIMENSION(:), POINTER :: temp
	temp=>sa%irow
	sa%irow=>sa%jcol
	sa%jcol=>temp
	END SUBROUTINE sprstp
