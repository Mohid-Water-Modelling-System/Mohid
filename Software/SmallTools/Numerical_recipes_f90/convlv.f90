	FUNCTION convlv(data,respns,isign)
	USE nrtype; USE nrutil, ONLY : assert,nrerror
	USE nr, ONLY : realft
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: data
	REAL(SP), DIMENSION(:), INTENT(IN) :: respns
	INTEGER(I4B), INTENT(IN) :: isign
	REAL(SP), DIMENSION(size(data)) :: convlv
	INTEGER(I4B) :: no2,n,m
	COMPLEX(SPC), DIMENSION(size(data)/2) :: tmpd,tmpr
	n=size(data)
	m=size(respns)
	call assert(iand(n,n-1)==0, 'n must be a power of 2 in convlv')
	call assert(mod(m,2)==1, 'm must be odd in convlv')
	convlv(1:m)=respns(:)
	convlv(n-(m-3)/2:n)=convlv((m+3)/2:m)
	convlv((m+3)/2:n-(m-1)/2)=0.0
	no2=n/2
	call realft(data,1,tmpd)
	call realft(convlv,1,tmpr)
	if (isign == 1) then
		tmpr(1)=cmplx(real(tmpd(1))*real(tmpr(1))/no2, &
			aimag(tmpd(1))*aimag(tmpr(1))/no2, kind=spc)
		tmpr(2:)=tmpd(2:)*tmpr(2:)/no2
	else if (isign == -1) then
		if (any(abs(tmpr(2:)) == 0.0) .or. real(tmpr(1)) == 0.0 &
			.or. aimag(tmpr(1)) == 0.0) call nrerror &
			('deconvolving at response zero in convlv')
		tmpr(1)=cmplx(real(tmpd(1))/real(tmpr(1))/no2, &
			aimag(tmpd(1))/aimag(tmpr(1))/no2, kind=spc)
		tmpr(2:)=tmpd(2:)/tmpr(2:)/no2
	else
		call nrerror('no meaning for isign in convlv')
	end if
	call realft(convlv,-1,tmpr)
	END FUNCTION convlv
