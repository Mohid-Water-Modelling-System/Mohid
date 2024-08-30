	FUNCTION erfc_s(x)
	USE nrtype
	USE nr, ONLY : gammp,gammq
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: erfc_s
	erfc_s=merge(1.0_sp+gammp(0.5_sp,x**2),gammq(0.5_sp,x**2), x < 0.0)
	END FUNCTION erfc_s


	FUNCTION erfc_v(x)
	USE nrtype
	USE nr, ONLY : gammp,gammq
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(size(x)) :: erfc_v
	LOGICAL(LGT), DIMENSION(size(x)) :: mask
	mask = (x < 0.0)
	erfc_v=merge(1.0_sp+gammp(spread(0.5_sp,1,size(x)), &
		merge(x,0.0_sp,mask)**2),gammq(spread(0.5_sp,1,size(x)), &
		merge(x,0.0_sp,.not. mask)**2),mask)
	END FUNCTION erfc_v
