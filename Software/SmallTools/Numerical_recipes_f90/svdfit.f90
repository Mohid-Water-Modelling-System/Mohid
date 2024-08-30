	SUBROUTINE svdfit(x,y,sig,a,v,w,chisq,funcs)
	USE nrtype; USE nrutil, ONLY : assert_eq,vabs
	USE nr, ONLY : svbksb,svdcmp
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x,y,sig
	REAL(SP), DIMENSION(:), INTENT(OUT) :: a,w
	REAL(SP), DIMENSION(:,:), INTENT(OUT) :: v
	REAL(SP), INTENT(OUT) :: chisq
	INTERFACE
		FUNCTION funcs(x,n)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		INTEGER(I4B), INTENT(IN) :: n
		REAL(SP), DIMENSION(n) :: funcs
		END FUNCTION funcs
	END INTERFACE
	REAL(SP), PARAMETER :: TOL=1.0e-5_sp
	INTEGER(I4B) :: i,ma,n
	REAL(SP), DIMENSION(size(x)) :: b,sigi
	REAL(SP), DIMENSION(size(x),size(a)) :: u,usav
	n=assert_eq(size(x),size(y),size(sig),'svdfit: n')
	ma=assert_eq(size(a),size(v,1),size(v,2),size(w),'svdfit: ma')
	sigi=1.0_sp/sig
	b=y*sigi
	do i=1,n
		usav(i,:)=funcs(x(i),ma)
	end do
	u=usav*spread(sigi,dim=2,ncopies=ma)
	usav=u
	call svdcmp(u,w,v)
	where (w < TOL*maxval(w)) w=0.0
	call svbksb(u,w,v,b,a)
	chisq=vabs(matmul(usav,a)-b)**2
	END SUBROUTINE svdfit
