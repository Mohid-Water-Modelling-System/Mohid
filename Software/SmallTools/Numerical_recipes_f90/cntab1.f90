	SUBROUTINE cntab1(nn,chisq,df,prob,cramrv,ccc)
	USE nrtype; USE nrutil, ONLY : outerprod
	USE nr, ONLY : gammq
	IMPLICIT NONE
	INTEGER(I4B), DIMENSION(:,:), INTENT(IN) :: nn
	REAL(SP), INTENT(OUT) :: chisq,df,prob,cramrv,ccc
	REAL(SP), PARAMETER :: TINY=1.0e-30_sp
	INTEGER(I4B) :: nni,nnj
	REAL(SP) :: sumn
	REAL(SP), DIMENSION(size(nn,1)) :: sumi
	REAL(SP), DIMENSION(size(nn,2)) :: sumj
	REAL(SP), DIMENSION(size(nn,1),size(nn,2)) :: expctd
	sumi(:)=sum(nn(:,:),dim=2)
	sumj(:)=sum(nn(:,:),dim=1)
	sumn=sum(sumi(:))
	nni=size(sumi)-count(sumi(:) == 0.0)
	nnj=size(sumj)-count(sumj(:) == 0.0)
	df=nni*nnj-nni-nnj+1
	expctd(:,:)=outerprod(sumi(:),sumj(:))/sumn
	chisq=sum((nn(:,:)-expctd(:,:))**2/(expctd(:,:)+TINY))
	prob=gammq(0.5_sp*df,0.5_sp*chisq)
	cramrv=sqrt(chisq/(sumn*min(nni-1,nnj-1)))
	ccc=sqrt(chisq/(chisq+sumn))
	END SUBROUTINE cntab1
