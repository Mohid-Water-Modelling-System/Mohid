	FUNCTION savgol(nl,nrr,ld,m)
	USE nrtype; USE nrutil, ONLY : arth,assert,poly
	USE nr, ONLY : lubksb,ludcmp
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: nl,nrr,ld,m
	REAL(SP), DIMENSION(nl+nrr+1) :: savgol
	INTEGER(I4B) :: imj,ipj,mm,np
	INTEGER(I4B), DIMENSION(m+1) :: indx
	REAL(SP) :: d,sm
	REAL(SP), DIMENSION(m+1) :: b
	REAL(SP), DIMENSION(m+1,m+1) :: a
	INTEGER(I4B) :: irng(nl+nrr+1)
	call assert(nl >= 0, nrr >= 0, ld <= m, nl+nrr >= m, 'savgol args')
	do ipj=0,2*m
		sm=sum(arth(1.0_sp,1.0_sp,nrr)**ipj)+&
			sum(arth(-1.0_sp,-1.0_sp,nl)**ipj)
		if (ipj == 0) sm=sm+1.0_sp
		mm=min(ipj,2*m-ipj)
		do imj=-mm,mm,2
			a(1+(ipj+imj)/2,1+(ipj-imj)/2)=sm
		end do
	end do
	call ludcmp(a(:,:),indx(:),d)
	b(:)=0.0
	b(ld+1)=1.0
	call lubksb(a(:,:),indx(:),b(:))
	savgol(:)=0.0
	irng(:)=arth(-nl,1,nrr+nl+1)
	np=nl+nrr+1
	savgol(mod(np-irng(:),np)+1)=poly(real(irng(:),sp),b(:))
	END FUNCTION savgol
