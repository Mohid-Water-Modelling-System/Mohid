	SUBROUTINE orthog(anu,alpha,beta,a,b)
	USE nrtype; USE nrutil, ONLY : assert_eq
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: anu,alpha,beta
	REAL(SP), DIMENSION(:), INTENT(OUT) :: a,b
	INTEGER(I4B) :: k,n,ndum
	REAL(SP), DIMENSION(2*size(a)+1,2*size(a)+1) :: sig
	n=assert_eq(size(a),size(b),'orthog: n')
	ndum=assert_eq(2*n,size(alpha)+1,size(anu),size(beta)+1,'orthog: ndum')
	sig(1,3:2*n)=0.0
	sig(2,2:2*n+1)=anu(1:2*n)
	a(1)=alpha(1)+anu(2)/anu(1)
	b(1)=0.0
	do k=3,n+1
		sig(k,k:2*n-k+3)=sig(k-1,k+1:2*n-k+4)+(alpha(k-1:2*n-k+2) &
			-a(k-2))*sig(k-1,k:2*n-k+3)-b(k-2)*sig(k-2,k:2*n-k+3) &
			+beta(k-1:2*n-k+2)*sig(k-1,k-1:2*n-k+2)
		a(k-1)=alpha(k-1)+sig(k,k+1)/sig(k,k)-sig(k-1,k)/sig(k-1,k-1)
		b(k-1)=sig(k,k)/sig(k-1,k-1)
	end do
	END SUBROUTINE orthog
