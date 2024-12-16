	SUBROUTINE mpinv(u,v,n,m)
	USE nrtype; USE nrutil, ONLY : poly
	USE nr, ONLY : mpmul
	USE mpops, ONLY : mpmov,mpneg
	IMPLICIT NONE
	CHARACTER(1), DIMENSION(:), INTENT(OUT) :: u
	CHARACTER(1), DIMENSION(:), INTENT(IN) :: v
	INTEGER(I4B), INTENT(IN) :: n,m
	INTEGER(I4B), PARAMETER :: MF=4
	REAL(SP), PARAMETER :: BI=1.0_sp/256.0_sp
	INTEGER(I4B) :: i,j,mm
	REAL(SP) :: fu
	CHARACTER(1), DIMENSION(:), ALLOCATABLE :: rr,s
	allocate(rr(max(n,m)+n+1),s(n))
	mm=min(MF,m)
	fu=1.0_sp/poly(BI,real(ichar(v(:)),sp))
	do j=1,n
		i=int(fu)
		u(j)=char(i)
		fu=256.0_sp*(fu-i)
	end do
	do
		call mpmul(rr,u,v,n,m)
		call mpmov(s,rr(2:),n)
		call mpneg(s,n)
		s(1)=char(ichar(s(1))-254)
		call mpmul(rr,s,u,n,n)
		call mpmov(u,rr(2:),n)
		if (all(ichar(s(2:n-1)) == 0)) exit
	end do
	deallocate(rr,s)
	END SUBROUTINE mpinv
