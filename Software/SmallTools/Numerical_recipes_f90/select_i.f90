	FUNCTION select_inplace(k,arr)
	USE nrtype
	USE nr, ONLY : select
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: k
	REAL(SP), DIMENSION(:), INTENT(IN) :: arr
	REAL(SP) :: select_inplace
	REAL(SP), DIMENSION(size(arr)) :: tarr
	tarr=arr
	select_inplace=select(k,tarr)
	END FUNCTION select_inplace
