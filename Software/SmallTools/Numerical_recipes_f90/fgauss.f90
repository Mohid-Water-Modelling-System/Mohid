	SUBROUTINE fgauss(x,a,y,dyda)
	USE nrtype; USE nrutil, ONLY : assert_eq
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x,a
	REAL(SP), DIMENSION(:), INTENT(OUT) :: y
	REAL(SP), DIMENSION(:,:), INTENT(OUT) :: dyda
	INTEGER(I4B) :: i,na,nx
	REAL(SP), DIMENSION(size(x)) :: arg,ex,fac
	nx=assert_eq(size(x),size(y),size(dyda,1),'fgauss: nx')
	na=assert_eq(size(a),size(dyda,2),'fgauss: na')
	y(:)=0.0
	do i=1,na-1,3
		arg(:)=(x(:)-a(i+1))/a(i+2)
		ex(:)=exp(-arg(:)**2)
		fac(:)=a(i)*ex(:)*2.0_sp*arg(:)
		y(:)=y(:)+a(i)*ex(:)
		dyda(:,i)=ex(:)
		dyda(:,i+1)=fac(:)/a(i+2)
		dyda(:,i+2)=fac(:)*arg(:)/a(i+2)
	end do
	END SUBROUTINE fgauss
