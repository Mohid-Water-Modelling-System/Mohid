	FUNCTION interp(uc)
	USE nrtype; USE nrutil, ONLY : assert_eq
	IMPLICIT NONE
	REAL(DP), DIMENSION(:,:), INTENT(IN) :: uc
	REAL(DP), DIMENSION(2*size(uc,1)-1,2*size(uc,1)-1) :: interp
	INTEGER(I4B) :: nc,nf
	nc=assert_eq(size(uc,1),size(uc,2),'interp')
	nf=2*nc-1
	interp(1:nf:2,1:nf:2)=uc(1:nc,1:nc)
	interp(2:nf-1:2,1:nf:2)=0.5_dp*(interp(3:nf:2,1:nf:2)+ &
		interp(1:nf-2:2,1:nf:2))
	interp(1:nf,2:nf-1:2)=0.5_dp*(interp(1:nf,3:nf:2)+interp(1:nf,1:nf-2:2))
	END FUNCTION interp
