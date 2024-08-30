	FUNCTION ratval_s(x,cof,mm,kk)
	USE nrtype; USE nrutil, ONLY : poly
	IMPLICIT NONE
	REAL(DP), INTENT(IN) :: x
	INTEGER(I4B), INTENT(IN) :: mm,kk
	REAL(DP), DIMENSION(mm+kk+1), INTENT(IN) :: cof
	REAL(DP) :: ratval_s
	ratval_s=poly(x,cof(1:mm+1))/(1.0_dp+x*poly(x,cof(mm+2:mm+kk+1)))
	END FUNCTION ratval_s

	FUNCTION ratval_v(x,cof,mm,kk)
	USE nrtype; USE nrutil, ONLY : poly
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), INTENT(IN) :: x
	INTEGER(I4B), INTENT(IN) :: mm,kk
	REAL(DP), DIMENSION(mm+kk+1), INTENT(IN) :: cof
	REAL(DP), DIMENSION(size(x)) :: ratval_v
	ratval_v=poly(x,cof(1:mm+1))/(1.0_dp+x*poly(x,cof(mm+2:mm+kk+1)))
	END FUNCTION ratval_v
