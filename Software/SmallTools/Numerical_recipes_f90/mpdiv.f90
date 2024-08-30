	SUBROUTINE mpdiv(q,r,u,v,n,m)
	USE nrtype; USE nrutil, ONLY : nrerror
	USE nr, ONLY : mpinv,mpmul
	USE mpops, ONLY : mpsad,mpmov,mpsub
	IMPLICIT NONE
	CHARACTER(1), DIMENSION(:), INTENT(OUT) :: q,r
	CHARACTER(1), DIMENSION(:), INTENT(IN) :: u,v
!	The logical dimensions are: CHARACTER(1) :: q(n-m+1),r(m),u(n),v(m)
	INTEGER(I4B), INTENT(IN) :: n,m
	INTEGER(I4B), PARAMETER :: MACC=6
	INTEGER(I4B) :: is
	CHARACTER(1), DIMENSION(:), ALLOCATABLE, TARGET :: rr,s
	CHARACTER(1), DIMENSION(:), POINTER :: rr2,s3
	allocate(rr(2*(n+MACC)),s(2*(n+MACC)))
	rr2=>rr(2:)
	s3=>s(3:)
	call mpinv(s,v,n+MACC,m)
	call mpmul(rr,s,u,n+MACC,n)
	call mpsad(s,rr,n+n+MACC/2,1)
	call mpmov(q,s3,n-m+1)
	call mpmul(rr,q,v,n-m+1,m)
	call mpsub(is,rr2,u,rr2,n)
	if (is /= 0) call nrerror('MACC too small in mpdiv')
	call mpmov(r,rr(n-m+2:),m)
	deallocate(rr,s)
	END SUBROUTINE mpdiv
