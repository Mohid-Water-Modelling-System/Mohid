	SUBROUTINE ksone(data,func,d,prob)
	USE nrtype; USE nrutil, ONLY : arth
	USE nr, ONLY : probks,sort
	IMPLICIT NONE
	REAL(SP), INTENT(OUT) :: d,prob
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: data
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER(I4B) :: n
	REAL(SP) :: en
	REAL(SP), DIMENSION(size(data)) :: fvals
	REAL(SP), DIMENSION(size(data)+1) :: temp
	call sort(data)
	n=size(data)
	en=n
	fvals(:)=func(data(:))
	temp=arth(0,1,n+1)/en
	d=maxval(max(abs(temp(1:n)-fvals(:)), &
		abs(temp(2:n+1)-fvals(:))))
	en=sqrt(en)
	prob=probks((en+0.12_sp+0.11_sp/en)*d)
	END SUBROUTINE ksone
