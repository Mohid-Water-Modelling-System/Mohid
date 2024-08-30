	SUBROUTINE gaujac(x,w,alf,bet)
	USE nrtype; USE nrutil, ONLY : arth,assert_eq,nrerror
	USE nr, ONLY : gammln
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: alf,bet
	REAL(SP), DIMENSION(:), INTENT(OUT) :: x,w
	REAL(DP), PARAMETER :: EPS=3.0e-14_dp
	INTEGER(I4B) :: its,j,n
	INTEGER(I4B), PARAMETER :: MAXIT=10
	REAL(DP) :: alfbet,a,c,temp
	REAL(DP), DIMENSION(size(x)) :: b,p1,p2,p3,pp,z,z1
	LOGICAL(LGT), DIMENSION(size(x)) :: unfinished
	n=assert_eq(size(x),size(w),'gaujac')
	alfbet=alf+bet
	z=cos(PI*(arth(1,1,n)-0.25_dp+0.5_dp*alf)/(n+0.5_dp*(alfbet+1.0_dp)))
	unfinished=.true.
	do its=1,MAXIT
		temp=2.0_dp+alfbet
		where (unfinished)
			p1=(alf-bet+temp*z)/2.0_dp
			p2=1.0
		end where
		do j=2,n
			a=2*j*(j+alfbet)*temp
			temp=temp+2.0_dp
			c=2.0_dp*(j-1.0_dp+alf)*(j-1.0_dp+bet)*temp
			where (unfinished)
				p3=p2
				p2=p1
				b=(temp-1.0_dp)*(alf*alf-bet*bet+temp*&
					(temp-2.0_dp)*z)
				p1=(b*p2-c*p3)/a
			end where
		end do
		where (unfinished)
			pp=(n*(alf-bet-temp*z)*p1+2.0_dp*(n+alf)*&
				(n+bet)*p2)/(temp*(1.0_dp-z*z))
			z1=z
			z=z1-p1/pp
			unfinished=(abs(z-z1) > EPS)
		end where
		if (.not. any(unfinished)) exit
	end do
	if (its == MAXIT+1) call nrerror('too many iterations in gaujac')
	x=z
	w=exp(gammln(alf+n)+gammln(bet+n)-gammln(n+1.0_sp)-&
		gammln(n+alf+bet+1.0_sp))*temp*2.0_sp**alfbet/(pp*p2)
	END SUBROUTINE gaujac
