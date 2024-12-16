	SUBROUTINE qrdcmp(a,c,d,sing)
	USE nrtype; USE nrutil, ONLY : assert_eq,outerprod,vabs
	IMPLICIT NONE
	REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
	REAL(SP), DIMENSION(:), INTENT(OUT) :: c,d
	LOGICAL(LGT), INTENT(OUT) :: sing
	INTEGER(I4B) :: k,n
	REAL(SP) :: scale,sigma
	n=assert_eq(size(a,1),size(a,2),size(c),size(d),'qrdcmp')
	sing=.false.
	do k=1,n-1
		scale=maxval(abs(a(k:n,k)))
		if (scale == 0.0) then
			sing=.true.
			c(k)=0.0
			d(k)=0.0
		else
			a(k:n,k)=a(k:n,k)/scale
			sigma=sign(vabs(a(k:n,k)),a(k,k))
			a(k,k)=a(k,k)+sigma
			c(k)=sigma*a(k,k)
			d(k)=-scale*sigma
			a(k:n,k+1:n)=a(k:n,k+1:n)-outerprod(a(k:n,k),&
				matmul(a(k:n,k),a(k:n,k+1:n)))/c(k)
		end if
	end do
	d(n)=a(n,n)
	if (d(n) == 0.0) sing=.true.
	END SUBROUTINE qrdcmp
