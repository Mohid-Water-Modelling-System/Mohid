	SUBROUTINE pwt(a,isign)
	USE nrtype; USE nrutil, ONLY : arth,nrerror
	USE pwtcom
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: a
	INTEGER(I4B), INTENT(IN) :: isign
	REAL(SP), DIMENSION(size(a)) :: wksp
	INTEGER(I4B), DIMENSION(size(a)/2) :: jf,jr
	INTEGER(I4B) :: k,n,nh,nmod
	n=size(a)
	if (n < 4) RETURN
	if (ncof == 0) call nrerror('pwt: must call pwtset before pwt')
	nmod=ncof*n
	nh=n/2
	wksp(:)=0.0
	jf=iand(n-1,arth(2+nmod+ioff,2,nh))
	jr=iand(n-1,arth(2+nmod+joff,2,nh))
	do k=1,ncof
		if (isign >= 0) then
			wksp(1:nh)=wksp(1:nh)+cc(k)*a(jf+1)
			wksp(nh+1:n)=wksp(nh+1:n)+cr(k)*a(jr+1)
		else
			wksp(jf+1)=wksp(jf+1)+cc(k)*a(1:nh)
			wksp(jr+1)=wksp(jr+1)+cr(k)*a(nh+1:n)
		end if
		if (k == ncof) exit
		jf=iand(n-1,jf+1)
		jr=iand(n-1,jr+1)
	end do
	a(:)=wksp(:)
	END SUBROUTINE pwt
