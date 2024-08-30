	SUBROUTINE sort3(arr,slave1,slave2)
	USE nrtype; USE nrutil, ONLY : assert_eq
	USE nr, ONLY : indexx
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr,slave1,slave2
	INTEGER(I4B) :: ndum
	INTEGER(I4B), DIMENSION(size(arr)) :: index
	ndum=assert_eq(size(arr),size(slave1),size(slave2),'sort3')
	call indexx(arr,index)
	arr=arr(index)
	slave1=slave1(index)
	slave2=slave2(index)
	END SUBROUTINE sort3
