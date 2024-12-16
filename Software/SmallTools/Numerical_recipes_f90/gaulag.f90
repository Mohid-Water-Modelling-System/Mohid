	SUBROUTINE gaulag(x,w,alf)
	USE nrtype; USE nrutil, ONLY : arth,assert_eq,nrerror
	USE nr, ONLY : gammln
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: alf
	REAL(SP), DIMENSION(:), INTENT(OUT) :: x,w
	REAL(DP), PARAMETER :: EPS=3.0e-13_dp
	INTEGER(I4B) :: its,j,n
	INTEGER(I4B), PARAMETER :: MAXIT=10
	REAL(SP) :: anu
	REAL(SP), PARAMETER :: C1=9.084064e-01_sp,C2=5.214976e-02_sp,&
		C3=2.579930e-03_sp,C4=3.986126e-03_sp
	REAL(SP), DIMENSION(size(x)) :: rhs,r2,r3,theta
	REAL(DP), DIMENSION(size(x)) :: p1,p2,p3,pp,z,z1
	LOGICAL(LGT), DIMENSION(size(x)) :: unfinished
	n=assert_eq(size(x),size(w),'gaulag')
	anu=4.0_sp*n+2.0_sp*alf+2.0_sp
	rhs=arth(4*n-1,-4,n)*PI/anu
	r3=rhs**(1.0_sp/3.0_sp)
	r2=r3**2
	theta=r3*(C1+r2*(C2+r2*(C3+r2*C4)))
	z=anu*cos(theta)**2
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
				p1=((2.0_dp*j-1.0_dp+alf-z)*p2-(j-1.0_dp+alf)*p3)/j
			end where
		end do
		where (unfinished)
			pp=(n*p1-(n+alf)*p2)/z
			z1=z
			z=z1-p1/pp
			unfinished=(abs(z-z1) > EPS*z)
		end where
		if (.not. any(unfinished)) exit
	end do
	if (its == MAXIT+1) call nrerror('too many iterations in gaulag')
	x=z
	w=-exp(gammln(alf+n)-gammln(real(n,sp)))/(pp*n*p2)
	END SUBROUTINE gaulag
