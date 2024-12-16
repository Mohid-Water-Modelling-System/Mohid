	SUBROUTINE fitexy(x,y,sigx,sigy,a,b,siga,sigb,chi2,q)
	USE nrtype; USE nrutil, ONLY : assert_eq,swap
	USE nr, ONLY : avevar,brent,fit,gammq,mnbrak,zbrent
	USE chixyfit
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x,y,sigx,sigy
	REAL(SP), INTENT(OUT) :: a,b,siga,sigb,chi2,q
	REAL(SP), PARAMETER :: POTN=1.571000_sp,BIG=1.0e30_sp,ACC=1.0e-3_sp
	INTEGER(I4B) :: j,n
	REAL(SP), DIMENSION(size(x)), TARGET :: xx,yy,sx,sy,ww
	REAL(SP), DIMENSION(6) :: ang,ch
	REAL(SP) :: amx,amn,varx,vary,scale,bmn,bmx,d1,d2,r2,&
		dum1,dum2,dum3,dum4,dum5
	n=assert_eq(size(x),size(y),size(sigx),size(sigy),'fitexy')
	xxp=>xx
	yyp=>yy
	sxp=>sx
	syp=>sy
	wwp=>ww
	call avevar(x,dum1,varx)
	call avevar(y,dum1,vary)
	scale=sqrt(varx/vary)
	xx(:)=x(:)
	yy(:)=y(:)*scale
	sx(:)=sigx(:)
	sy(:)=sigy(:)*scale
	ww(:)=sqrt(sx(:)**2+sy(:)**2)
	call fit(xx,yy,dum1,b,dum2,dum3,dum4,dum5,ww)
	offs=0.0
	ang(1)=0.0
	ang(2)=atan(b)
	ang(4)=0.0
	ang(5)=ang(2)
	ang(6)=POTN
	do j=4,6
		ch(j)=chixy(ang(j))
	end do
	call mnbrak(ang(1),ang(2),ang(3),ch(1),ch(2),ch(3),chixy)
	chi2=brent(ang(1),ang(2),ang(3),chixy,ACC,b)
	chi2=chixy(b)
	a=aa
	q=gammq(0.5_sp*(n-2),0.5_sp*chi2)
	r2=1.0_sp/sum(ww(:))
	bmx=BIG
	bmn=BIG
	offs=chi2+1.0_sp
	do j=1,6
		if (ch(j) > offs) then
			d1=mod(abs(ang(j)-b),PI)
			d2=PI-d1
			if (ang(j) < b) call swap(d1,d2)
			if (d1 < bmx) bmx=d1
			if (d2 < bmn) bmn=d2
		end if
	end do
	if (bmx <  BIG) then
		bmx=zbrent(chixy,b,b+bmx,ACC)-b
		amx=aa-a
		bmn=zbrent(chixy,b,b-bmn,ACC)-b
		amn=aa-a
		sigb=sqrt(0.5_sp*(bmx**2+bmn**2))/(scale*cos(b)**2)
		siga=sqrt(0.5_sp*(amx**2+amn**2)+r2)/scale
	else
		sigb=BIG
		siga=BIG
	end if
	a=a/scale
	b=tan(b)/scale
	END SUBROUTINE fitexy
