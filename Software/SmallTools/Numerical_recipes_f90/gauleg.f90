	SUBROUTINE gauleg(x1,x2,x,w)
	USE nrtype; USE nrutil, ONLY : arth,assert_eq,nrerror
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x1,x2
	REAL(SP), DIMENSION(:), INTENT(OUT) :: x,w
	REAL(DP), PARAMETER :: EPS=3.0e-14_dp
	INTEGER(I4B) :: its,j,m,n
	INTEGER(I4B), PARAMETER :: MAXIT=10
	REAL(DP) :: xl,xm
	REAL(DP), DIMENSION((size(x)+1)/2) :: p1,p2,p3,pp,z,z1
	LOGICAL(LGT), DIMENSION((size(x)+1)/2) :: unfinished
	n=assert_eq(size(x),size(w),'gauleg')
	m=(n+1)/2
	xm=0.5_dp*(x2+x1)
	xl=0.5_dp*(x2-x1)
	z=cos(PI_D*(arth(1,1,m)-0.25_dp)/(n+0.5_dp))
	unfinished=.true.
	do its=1,MAXIT
		where (unfinished)
			p1=1.0
			p2=0.0
		end where
		do j=1,n
			where (unfinished)
				p3=p2
				p2=p1
				p1=((2.0_dp*j-1.0_dp)*z*p2-(j-1.0_dp)*p3)/j
			end where
		end do
		where (unfinished)
			pp=n*(z*p1-p2)/(z*z-1.0_dp)
			z1=z
			z=z1-p1/pp
			unfinished=(abs(z-z1) > EPS)
		end where
		if (.not. any(unfinished)) exit
	end do
	if (its == MAXIT+1) call nrerror('too many iterations in gauleg')
	x(1:m)=xm-xl*z
	x(n:n-m+1:-1)=xm+xl*z
	w(1:m)=2.0_dp*xl/((1.0_dp-z**2)*pp**2)
	w(n:n-m+1:-1)=w(1:m)
	END SUBROUTINE gauleg
