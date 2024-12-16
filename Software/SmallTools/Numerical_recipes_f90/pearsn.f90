	SUBROUTINE pearsn(x,y,r,prob,z)
	USE nrtype; USE nrutil, ONLY : assert_eq
	USE nr, ONLY : betai
	IMPLICIT NONE
	REAL(SP), INTENT(OUT) :: r,prob,z
	REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
	REAL(SP), PARAMETER :: TINY=1.0e-20_sp
	REAL(SP), DIMENSION(size(x)) :: xt,yt
	REAL(SP) :: ax,ay,df,sxx,sxy,syy,t
	INTEGER(I4B) :: n
	n=assert_eq(size(x),size(y),'pearsn')
	ax=sum(x)/n
	ay=sum(y)/n
	xt(:)=x(:)-ax
	yt(:)=y(:)-ay
	sxx=dot_product(xt,xt)
	syy=dot_product(yt,yt)
	sxy=dot_product(xt,yt)
	r=sxy/(sqrt(sxx*syy)+TINY)
	z=0.5_sp*log(((1.0_sp+r)+TINY)/((1.0_sp-r)+TINY))
	df=n-2
	t=r*sqrt(df/(((1.0_sp-r)+TINY)*((1.0_sp+r)+TINY)))
	prob=betai(0.5_sp*df,0.5_sp,df/(df+t**2))
!	prob=erfcc(abs(z*sqrt(n-1.0_sp))/SQRT2)
	END SUBROUTINE pearsn
