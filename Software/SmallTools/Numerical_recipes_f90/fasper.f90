	SUBROUTINE fasper(x,y,ofac,hifac,px,py,jmax,prob)
	USE nrtype; USE nrutil, ONLY : arth,assert_eq,imaxloc,nrerror
	USE nr, ONLY : avevar,realft
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
	REAL(SP), INTENT(IN) :: ofac,hifac
	INTEGER(I4B), INTENT(OUT) :: jmax
	REAL(SP), INTENT(OUT) :: prob
	REAL(SP), DIMENSION(:), POINTER :: px,py
	INTEGER(I4B), PARAMETER :: MACC=4
	INTEGER(I4B) :: j,k,n,ndim,nfreq,nfreqt,nout
	REAL(SP) :: ave,ck,ckk,cterm,cwt,den,df,effm,expy,fac,fndim,hc2wt,&
		hs2wt,hypo,sterm,swt,var,xdif,xmax,xmin
	REAL(SP), DIMENSION(:), ALLOCATABLE :: wk1,wk2
	LOGICAL(LGT), SAVE :: init=.true.
	n=assert_eq(size(x),size(y),'fasper')
	if (init) then
		init=.false.
		nullify(px,py)
	else
		if (associated(px)) deallocate(px)
		if (associated(py)) deallocate(py)
	end if
	nfreqt=ofac*hifac*n*MACC
	nfreq=64
	do
		if (nfreq >= nfreqt) exit
		nfreq=nfreq*2
	end do
	ndim=2*nfreq
	allocate(wk1(ndim),wk2(ndim))
	call avevar(y(1:n),ave,var)
	xmax=maxval(x(:))
	xmin=minval(x(:))
	xdif=xmax-xmin
	wk1(1:ndim)=0.0
	wk2(1:ndim)=0.0
	fac=ndim/(xdif*ofac)
	fndim=ndim
	do j=1,n
		ck=1.0_sp+mod((x(j)-xmin)*fac,fndim)
		ckk=1.0_sp+mod(2.0_sp*(ck-1.0_sp),fndim)
		call spreadval(y(j)-ave,wk1,ck,MACC)
		call spreadval(1.0_sp,wk2,ckk,MACC)
	end do
	call realft(wk1(1:ndim),1)
	call realft(wk2(1:ndim),1)
	df=1.0_sp/(xdif*ofac)
	nout=0.5_sp*ofac*hifac*n
	allocate(px(nout),py(nout))
	k=3
	do j=1,nout
		hypo=sqrt(wk2(k)**2+wk2(k+1)**2)
		hc2wt=0.5_sp*wk2(k)/hypo
		hs2wt=0.5_sp*wk2(k+1)/hypo
		cwt=sqrt(0.5_sp+hc2wt)
		swt=sign(sqrt(0.5_sp-hc2wt),hs2wt)
		den=0.5_sp*n+hc2wt*wk2(k)+hs2wt*wk2(k+1)
		cterm=(cwt*wk1(k)+swt*wk1(k+1))**2/den
		sterm=(cwt*wk1(k+1)-swt*wk1(k))**2/(n-den)
		px(j)=j*df
		py(j)=(cterm+sterm)/(2.0_sp*var)
		k=k+2
	end do
	deallocate(wk1,wk2)
	jmax=imaxloc(py(1:nout))
	expy=exp(-py(jmax))
	effm=2.0_sp*nout/ofac
	prob=effm*expy
	if (prob > 0.01_sp) prob=1.0_sp-(1.0_sp-expy)**effm
	CONTAINS
!BL
	SUBROUTINE spreadval(y,yy,x,m)
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: y,x
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: yy
	INTEGER(I4B), INTENT(IN) :: m
	INTEGER(I4B) :: ihi,ilo,ix,j,nden,n
	REAL(SP) :: fac
	INTEGER(I4B), DIMENSION(10) :: nfac = (/ &
		1,1,2,6,24,120,720,5040,40320,362880 /)
	if (m > 10) call nrerror('factorial table too small in spreadval')
	n=size(yy)
	ix=x
	if (x == real(ix,sp)) then
		yy(ix)=yy(ix)+y
	else
		ilo=min(max(int(x-0.5_sp*m+1.0_sp),1),n-m+1)
		ihi=ilo+m-1
		nden=nfac(m)
		fac=product(x-arth(ilo,1,m))
		yy(ihi)=yy(ihi)+y*fac/(nden*(x-ihi))
		do j=ihi-1,ilo,-1
			nden=(nden/(j+1-ilo))*(j-ihi)
			yy(j)=yy(j)+y*fac/(nden*(x-j))
		end do
	end if
	END SUBROUTINE spreadval
	END SUBROUTINE fasper
