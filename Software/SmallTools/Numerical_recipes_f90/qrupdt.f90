	SUBROUTINE qrupdt(r,qt,u,v)
	USE nrtype; USE nrutil, ONLY : assert_eq,ifirstloc
	USE nr, ONLY : rotate,pythag
	IMPLICIT NONE
	REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: r,qt
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: u
	REAL(SP), DIMENSION(:), INTENT(IN) :: v
	INTEGER(I4B) :: i,k,n
	n=assert_eq((/size(r,1),size(r,2),size(qt,1),size(qt,2),size(u),&
		size(v)/),'qrupdt')
	k=n+1-ifirstloc(u(n:1:-1) /= 0.0)
	if (k < 1) k=1
	do i=k-1,1,-1
		call rotate(r,qt,i,u(i),-u(i+1))
		u(i)=pythag(u(i),u(i+1))
	end do
	r(1,:)=r(1,:)+u(1)*v
	do i=1,k-1
		call rotate(r,qt,i,r(i,i),-r(i+1,i))
	end do
	END SUBROUTINE qrupdt
