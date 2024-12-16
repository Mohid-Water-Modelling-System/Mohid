	FUNCTION rank(index)
	USE nrtype; USE nrutil, ONLY : arth
	IMPLICIT NONE
	INTEGER(I4B), DIMENSION(:), INTENT(IN) :: index
	INTEGER(I4B), DIMENSION(size(index)) :: rank
	rank(index(:))=arth(1,1,size(index))
	END FUNCTION rank
