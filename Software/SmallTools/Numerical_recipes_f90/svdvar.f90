	SUBROUTINE svdvar(v,w,cvm)
	USE nrtype; USE nrutil, ONLY : assert_eq
	IMPLICIT NONE
	REAL(SP), DIMENSION(:,:), INTENT(IN) :: v
	REAL(SP), DIMENSION(:), INTENT(IN) :: w
	REAL(SP), DIMENSION(:,:), INTENT(OUT) :: cvm
	INTEGER(I4B) :: ma
	REAL(SP), DIMENSION(size(w)) :: wti
	ma=assert_eq((/size(v,1),size(v,2),size(w),size(cvm,1),size(cvm,2)/),&
		'svdvar')
	where (w /= 0.0)
		wti=1.0_sp/(w*w)
	elsewhere
		wti=0.0
	end where
	cvm=v*spread(wti,dim=1,ncopies=ma)
	cvm=matmul(cvm,transpose(v))
	END SUBROUTINE svdvar
