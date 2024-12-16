	SUBROUTINE moment(data,ave,adev,sdev,var,skew,curt)
	USE nrtype; USE nrutil, ONLY : nrerror
	IMPLICIT NONE
	REAL(SP), INTENT(OUT) :: ave,adev,sdev,var,skew,curt
	REAL(SP), DIMENSION(:), INTENT(IN) :: data
	INTEGER(I4B) :: n
	REAL(SP) :: ep
	REAL(SP), DIMENSION(size(data)) :: p,s
	n=size(data)
	if (n <= 1) call nrerror('moment: n must be at least 2')
	ave=sum(data(:))/n
	s(:)=data(:)-ave
	ep=sum(s(:))
	adev=sum(abs(s(:)))/n
	p(:)=s(:)*s(:)
	var=sum(p(:))
	p(:)=p(:)*s(:)
	skew=sum(p(:))
	p(:)=p(:)*s(:)
	curt=sum(p(:))
	var=(var-ep**2/n)/(n-1)
	sdev=sqrt(var)
	if (var /= 0.0) then
		skew=skew/(n*sdev**3)
		curt=curt/(n*var**2)-3.0_sp
	else
		call nrerror('moment: no skew or kurtosis when zero variance')
	end if
	END SUBROUTINE moment
