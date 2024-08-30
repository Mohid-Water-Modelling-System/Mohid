	FUNCTION rstrct(uf)
	USE nrtype; USE nrutil, ONLY : assert_eq
	IMPLICIT NONE
	REAL(DP), DIMENSION(:,:), INTENT(IN) :: uf
	REAL(DP), DIMENSION((size(uf,1)+1)/2,(size(uf,1)+1)/2) :: rstrct
	INTEGER(I4B) :: nc,nf
	nf=assert_eq(size(uf,1),size(uf,2),'rstrct')
	nc=(nf+1)/2
	rstrct(2:nc-1,2:nc-1)=0.5_dp*uf(3:nf-2:2,3:nf-2:2)+0.125_dp*(&
		uf(4:nf-1:2,3:nf-2:2)+uf(2:nf-3:2,3:nf-2:2)+&
		uf(3:nf-2:2,4:nf-1:2)+uf(3:nf-2:2,2:nf-3:2))
	rstrct(1:nc,1)=uf(1:nf:2,1)
	rstrct(1:nc,nc)=uf(1:nf:2,nf)
	rstrct(1,1:nc)=uf(1,1:nf:2)
	rstrct(nc,1:nc)=uf(nf,1:nf:2)
	END FUNCTION rstrct
