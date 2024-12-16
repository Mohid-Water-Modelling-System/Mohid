	FUNCTION toeplz(r,y)
	USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: r,y
	REAL(SP), DIMENSION(size(y)) :: toeplz
	INTEGER(I4B) :: m,m1,n,ndum
	REAL(SP) :: sd,sgd,sgn,shn,sxn
	REAL(SP), DIMENSION(size(y)) :: g,h,t
	n=size(y)
	ndum=assert_eq(2*n-1,size(r),'toeplz: ndum')
	if (r(n) == 0.0) call nrerror('toeplz: initial singular minor')
	toeplz(1)=y(1)/r(n)
	if (n == 1) RETURN
	g(1)=r(n-1)/r(n)
	h(1)=r(n+1)/r(n)
	do m=1,n
		m1=m+1
		sxn=-y(m1)+dot_product(r(n+1:n+m),toeplz(m:1:-1))
		sd=-r(n)+dot_product(r(n+1:n+m),g(1:m))
		if (sd == 0.0) exit
		toeplz(m1)=sxn/sd
		toeplz(1:m)=toeplz(1:m)-toeplz(m1)*g(m:1:-1)
		if (m1 == n) RETURN
		sgn=-r(n-m1)+dot_product(r(n-m:n-1),g(1:m))
		shn=-r(n+m1)+dot_product(r(n+m:n+1:-1),h(1:m))
		sgd=-r(n)+dot_product(r(n-m:n-1),h(m:1:-1))
		if (sd == 0.0 .or. sgd == 0.0) exit
		g(m1)=sgn/sgd
		h(m1)=shn/sd
		t(1:m)=g(1:m)
		g(1:m)=g(1:m)-g(m1)*h(m:1:-1)
		h(1:m)=h(1:m)-h(m1)*t(m:1:-1)
	end do
	if (m > n) call nrerror('toeplz: sanity check failed in routine')
	call nrerror('toeplz: singular principal minor')
	END FUNCTION toeplz
