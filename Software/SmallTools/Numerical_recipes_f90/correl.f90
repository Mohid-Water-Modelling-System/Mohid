	FUNCTION correl(data1,data2)
	USE nrtype; USE nrutil, ONLY : assert,assert_eq
	USE nr, ONLY : realft
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: data1,data2
	REAL(SP), DIMENSION(size(data1)) :: correl
	COMPLEX(SPC), DIMENSION(size(data1)/2) :: cdat1,cdat2
	INTEGER(I4B) :: no2,n
	n=assert_eq(size(data1),size(data2),'correl')
	call assert(iand(n,n-1)==0, 'n must be a power of 2 in correl')
	no2=n/2
	call realft(data1,1,cdat1)
	call realft(data2,1,cdat2)
	cdat1(1)=cmplx(real(cdat1(1))*real(cdat2(1))/no2, &
		aimag(cdat1(1))*aimag(cdat2(1))/no2, kind=spc)
	cdat1(2:)=cdat1(2:)*conjg(cdat2(2:))/no2
	call realft(correl,-1,cdat1)
	END FUNCTION correl
