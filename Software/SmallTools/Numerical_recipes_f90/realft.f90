	SUBROUTINE realft_sp(data,isign,zdata)
	USE nrtype; USE nrutil, ONLY : assert,assert_eq,zroots_unity
	USE nr, ONLY : four1
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: data
	INTEGER(I4B), INTENT(IN) :: isign
	COMPLEX(SPC), DIMENSION(:), OPTIONAL, TARGET :: zdata
	INTEGER(I4B) :: n,ndum,nh,nq
	COMPLEX(SPC), DIMENSION(size(data)/4) :: w
	COMPLEX(SPC), DIMENSION(size(data)/4-1) :: h1,h2
	COMPLEX(SPC), DIMENSION(:), POINTER :: cdata
	COMPLEX(SPC) :: z
	REAL(SP) :: c1=0.5_sp,c2
	n=size(data)
	call assert(iand(n,n-1)==0, 'n must be a power of 2 in realft_sp')
	nh=n/2
	nq=n/4
	if (present(zdata)) then
		ndum=assert_eq(n/2,size(zdata),'realft_sp')
		cdata=>zdata
		if (isign == 1) cdata=cmplx(data(1:n-1:2),data(2:n:2),kind=spc)
	else
		allocate(cdata(n/2))
		cdata=cmplx(data(1:n-1:2),data(2:n:2),kind=spc)
	end if
	if (isign == 1) then
		c2=-0.5_sp
		call four1(cdata,+1)
	else
		c2=0.5_sp
	end if
	w=zroots_unity(sign(n,isign),n/4)
	w=cmplx(-aimag(w),real(w),kind=spc)
	h1=c1*(cdata(2:nq)+conjg(cdata(nh:nq+2:-1)))
	h2=c2*(cdata(2:nq)-conjg(cdata(nh:nq+2:-1)))
	cdata(2:nq)=h1+w(2:nq)*h2
	cdata(nh:nq+2:-1)=conjg(h1-w(2:nq)*h2)
	z=cdata(1)
	if (isign == 1) then
		cdata(1)=cmplx(real(z)+aimag(z),real(z)-aimag(z),kind=spc)
	else
		cdata(1)=cmplx(c1*(real(z)+aimag(z)),c1*(real(z)-aimag(z)),kind=spc)
		call four1(cdata,-1)
	end if
	if (present(zdata)) then
		if (isign /= 1) then
			data(1:n-1:2)=real(cdata)
			data(2:n:2)=aimag(cdata)
		end if
	else
		data(1:n-1:2)=real(cdata)
		data(2:n:2)=aimag(cdata)
		deallocate(cdata)
	end if
	END SUBROUTINE realft_sp


	SUBROUTINE realft_dp(data,isign,zdata)
	USE nrtype; USE nrutil, ONLY : assert,assert_eq,zroots_unity
	USE nr, ONLY : four1
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), INTENT(INOUT) :: data
	INTEGER(I4B), INTENT(IN) :: isign
	COMPLEX(DPC), DIMENSION(:), OPTIONAL, TARGET :: zdata
	INTEGER(I4B) :: n,ndum,nh,nq
	COMPLEX(DPC), DIMENSION(size(data)/4) :: w
	COMPLEX(DPC), DIMENSION(size(data)/4-1) :: h1,h2
	COMPLEX(DPC), DIMENSION(:), POINTER :: cdata
	COMPLEX(DPC) :: z
	REAL(DP) :: c1=0.5_dp,c2
	n=size(data)
	call assert(iand(n,n-1)==0, 'n must be a power of 2 in realft_dp')
	nh=n/2
	nq=n/4
	if (present(zdata)) then
		ndum=assert_eq(n/2,size(zdata),'realft_dp')
		cdata=>zdata
		if (isign == 1) cdata=cmplx(data(1:n-1:2),data(2:n:2),kind=spc)
	else
		allocate(cdata(n/2))
		cdata=cmplx(data(1:n-1:2),data(2:n:2),kind=spc)
	end if
	if (isign == 1) then
		c2=-0.5_dp
		call four1(cdata,+1)
	else
		c2=0.5_dp
	end if
	w=zroots_unity(sign(n,isign),n/4)
	w=cmplx(-aimag(w),real(w),kind=dpc)
	h1=c1*(cdata(2:nq)+conjg(cdata(nh:nq+2:-1)))
	h2=c2*(cdata(2:nq)-conjg(cdata(nh:nq+2:-1)))
	cdata(2:nq)=h1+w(2:nq)*h2
	cdata(nh:nq+2:-1)=conjg(h1-w(2:nq)*h2)
	z=cdata(1)
	if (isign == 1) then
		cdata(1)=cmplx(real(z)+aimag(z),real(z)-aimag(z),kind=dpc)
	else
		cdata(1)=cmplx(c1*(real(z)+aimag(z)),c1*(real(z)-aimag(z)),kind=dpc)
		call four1(cdata,-1)
	end if
	if (present(zdata)) then
		if (isign /= 1) then
			data(1:n-1:2)=real(cdata)
			data(2:n:2)=aimag(cdata)
		end if
	else
		data(1:n-1:2)=real(cdata)
		data(2:n:2)=aimag(cdata)
		deallocate(cdata)
	end if
	END SUBROUTINE realft_dp

