	SUBROUTINE poldiv(u,v,q,r)
	USE nrtype; USE nrutil, ONLY : assert_eq
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: u,v
	REAL(SP), DIMENSION(:), INTENT(OUT) :: q,r
	INTEGER(I4B) :: i,n,nv
	n=assert_eq(size(u),size(q),size(r),'poldiv')
	nv=size(v)
	r(:)=u(:)
	q(:)=0.0
	do i=n-nv,0,-1
		q(i+1)=r(nv+i)/v(nv)
		r(i+1:nv+i-1)=r(i+1:nv+i-1)-q(i+1)*v(1:nv-1)
	end do
	r(nv:n)=0.0
	END SUBROUTINE poldiv
