	SUBROUTINE mpsqrt(w,u,v,n,m)
	USE nrtype; USE nrutil, ONLY : poly
	USE nr, ONLY : mpmul
	USE mpops, ONLY : mplsh,mpmov,mpneg,mpsdv
	IMPLICIT NONE
	CHARACTER(1), DIMENSION(:), INTENT(OUT) :: w,u
	CHARACTER(1), DIMENSION(:), INTENT(IN) :: v
	INTEGER(I4B), INTENT(IN) :: n,m
	INTEGER(I4B), PARAMETER :: MF=3
	REAL(SP), PARAMETER :: BI=1.0_sp/256.0_sp
	INTEGER(I4B) :: i,ir,j,mm
	REAL(SP) :: fu
	CHARACTER(1), DIMENSION(:), ALLOCATABLE :: r,s
	allocate(r(2*n),s(2*n))
	mm=min(m,MF)
	fu=1.0_sp/sqrt(poly(BI,real(ichar(v(:)),sp)))
	do j=1,n
		i=int(fu)
		u(j)=char(i)
		fu=256.0_sp*(fu-i)
	end do
	do
		call mpmul(r,u,u,n,n)
		call mplsh(r,n)
		call mpmul(s,r,v,n,min(m,n))
		call mplsh(s,n)
		call mpneg(s,n)
		s(1)=char(ichar(s(1))-253)
		call mpsdv(s,s,n,2,ir)
		if (any(ichar(s(2:n-1)) /= 0)) then
			call mpmul(r,s,u,n,n)
			call mpmov(u,r(2:),n)
			cycle
		end if
		call mpmul(r,u,v,n,min(m,n))
		call mpmov(w,r(2:),n)
		deallocate(r,s)
		RETURN
	end do
	END SUBROUTINE mpsqrt
