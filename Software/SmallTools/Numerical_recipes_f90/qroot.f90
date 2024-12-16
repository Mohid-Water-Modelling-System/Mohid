	SUBROUTINE qroot(p,b,c,eps)
	USE nrtype; USE nrutil, ONLY : nrerror
	USE nr, ONLY : poldiv
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: p
	REAL(SP), INTENT(INOUT) :: b,c
	REAL(SP), INTENT(IN) :: eps
	INTEGER(I4B), PARAMETER :: ITMAX=20
	REAL(SP), PARAMETER :: TINY=1.0e-6_sp
	INTEGER(I4B) :: iter,n
	REAL(SP) :: delb,delc,div,r,rb,rc,s,sb,sc
	REAL(SP), DIMENSION(3) :: d
	REAL(SP), DIMENSION(size(p)) :: q,qq,rem
	n=size(p)
	d(3)=1.0
	do iter=1,ITMAX
		d(2)=b
		d(1)=c
		call poldiv(p,d,q,rem)
		s=rem(1)
		r=rem(2)
		call poldiv(q(1:n-1),d(:),qq(1:n-1),rem(1:n-1))
		sc=-rem(1)
		rc=-rem(2)
		sb=-c*rc
		rb=sc-b*rc
		div=1.0_sp/(sb*rc-sc*rb)
		delb=(r*sc-s*rc)*div
		delc=(-r*sb+s*rb)*div
		b=b+delb
		c=c+delc
		if ((abs(delb) <= eps*abs(b) .or. abs(b) < TINY) .and. &
			(abs(delc) <= eps*abs(c) .or. abs(c) < TINY)) RETURN
	end do
	call nrerror('qroot: too many iterations')
	END SUBROUTINE qroot
