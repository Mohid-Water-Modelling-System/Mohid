	SUBROUTINE chsone(bins,ebins,knstrn,df,chsq,prob)
	USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
	USE nr, ONLY : gammq
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: knstrn
	REAL(SP), INTENT(OUT) :: df,chsq,prob
	REAL(SP), DIMENSION(:), INTENT(IN) :: bins,ebins
	INTEGER(I4B) :: ndum
	ndum=assert_eq(size(bins),size(ebins),'chsone')
	if (any(ebins(:) <= 0.0)) call nrerror('bad expected number in chsone')
	df=size(bins)-knstrn
	chsq=sum((bins(:)-ebins(:))**2/ebins(:))
	prob=gammq(0.5_sp*df,0.5_sp*chsq)
	END SUBROUTINE chsone
