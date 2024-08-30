	SUBROUTINE eulsum(sum,term,jterm)
	USE nrtype; USE nrutil, ONLY : poly_term,reallocate
	IMPLICIT NONE
	REAL(SP), INTENT(INOUT) :: sum
	REAL(SP), INTENT(IN) :: term
	INTEGER(I4B), INTENT(IN) :: jterm
	REAL(SP), DIMENSION(:), POINTER, SAVE :: wksp
	INTEGER(I4B), SAVE :: nterm
	LOGICAL(LGT), SAVE :: init=.true.
	if (init) then
		init=.false.
		nullify(wksp)
	end if
	if (jterm == 1) then
		nterm=1
		wksp=>reallocate(wksp,100)
		wksp(1)=term
		sum=0.5_sp*term
	else
		if (nterm+1 > size(wksp)) wksp=>reallocate(wksp,2*size(wksp))
		wksp(2:nterm+1)=0.5_sp*wksp(1:nterm)
		wksp(1)=term
		wksp(1:nterm+1)=poly_term(wksp(1:nterm+1),0.5_sp)
		if (abs(wksp(nterm+1)) <= abs(wksp(nterm))) then
			sum=sum+0.5_sp*wksp(nterm+1)
			nterm=nterm+1
		else
			sum=sum+wksp(nterm+1)
		end if
	end if
	END SUBROUTINE eulsum
