	SUBROUTINE cntab2(nn,h,hx,hy,hygx,hxgy,uygx,uxgy,uxy)
	USE nrtype
	IMPLICIT NONE
	INTEGER(I4B), DIMENSION(:,:), INTENT(IN) :: nn
	REAL(SP), INTENT(OUT) :: h,hx,hy,hygx,hxgy,uygx,uxgy,uxy
	REAL(SP), PARAMETER :: TINY=1.0e-30_sp
	REAL(SP) :: sumn
	REAL(SP), DIMENSION(size(nn,1)) :: sumi
	REAL(SP), DIMENSION(size(nn,2)) :: sumj
	sumi(:)=sum(nn(:,:),dim=2)
	sumj(:)=sum(nn(:,:),dim=1)
	sumn=sum(sumi(:))
	hx=-sum(sumi(:)*log(sumi(:)/sumn), mask=(sumi(:) /= 0.0) )/sumn
	hy=-sum(sumj(:)*log(sumj(:)/sumn), mask=(sumj(:) /= 0.0) )/sumn
	h=-sum(nn(:,:)*log(nn(:,:)/sumn), mask=(nn(:,:) /= 0) )/sumn
	hygx=h-hx
	hxgy=h-hy
	uygx=(hy-hygx)/(hy+TINY)
	uxgy=(hx-hxgy)/(hx+TINY)
	uxy=2.0_sp*(hx+hy-h)/(hx+hy+TINY)
	END SUBROUTINE cntab2
