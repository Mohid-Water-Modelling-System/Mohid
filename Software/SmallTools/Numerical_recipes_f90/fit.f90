	SUBROUTINE fit(x,y,a,b,siga,sigb,chi2,q,sig)
	USE nrtype; USE nrutil, ONLY : assert_eq
	USE nr, ONLY : gammq
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
	REAL(SP), INTENT(OUT) :: a,b,siga,sigb,chi2,q
	REAL(SP), DIMENSION(:), OPTIONAL, INTENT(IN) :: sig
	INTEGER(I4B) :: ndata
	REAL(SP) :: sigdat,ss,sx,sxoss,sy,st2
	REAL(SP), DIMENSION(size(x)), TARGET :: t
	REAL(SP), DIMENSION(:), POINTER :: wt
	if (present(sig)) then
		ndata=assert_eq(size(x),size(y),size(sig),'fit')
		wt=>t
		wt(:)=1.0_sp/(sig(:)**2)
		ss=sum(wt(:))
		sx=dot_product(wt,x)
		sy=dot_product(wt,y)
	else
		ndata=assert_eq(size(x),size(y),'fit')
		ss=real(size(x),sp)
		sx=sum(x)
		sy=sum(y)
	end if
	sxoss=sx/ss
	t(:)=x(:)-sxoss
	if (present(sig)) then
		t(:)=t(:)/sig(:)
		b=dot_product(t/sig,y)
	else
		b=dot_product(t,y)
	end if
	st2=dot_product(t,t)
	b=b/st2
	a=(sy-sx*b)/ss
	siga=sqrt((1.0_sp+sx*sx/(ss*st2))/ss)
	sigb=sqrt(1.0_sp/st2)
	t(:)=y(:)-a-b*x(:)
	q=1.0
	if (present(sig)) then
		t(:)=t(:)/sig(:)
		chi2=dot_product(t,t)
		if (ndata > 2) q=gammq(0.5_sp*(size(x)-2),0.5_sp*chi2)
	else
		chi2=dot_product(t,t)
		sigdat=sqrt(chi2/(size(x)-2))
		siga=siga*sigdat
		sigb=sigb*sigdat
	end if
	END SUBROUTINE fit
