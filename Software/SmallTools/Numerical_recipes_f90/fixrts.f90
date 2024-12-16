	SUBROUTINE fixrts(d)
	USE nrtype
	USE nr, ONLY : zroots
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: d
	INTEGER(I4B) :: i,m
	LOGICAL(LGT) :: polish
	COMPLEX(SPC), DIMENSION(size(d)+1) :: a
	COMPLEX(SPC), DIMENSION(size(d)) :: roots
	m=size(d)
	a(m+1)=cmplx(1.0_sp,kind=spc)
	a(m:1:-1)=cmplx(-d(1:m),kind=spc)
	polish=.true.
	call zroots(a(1:m+1),roots,polish)
	where (abs(roots) > 1.0) roots=1.0_sp/conjg(roots)
	a(1)=-roots(1)
	a(2:m+1)=cmplx(1.0_sp,kind=spc)
	do i=2,m
		a(2:i)=a(1:i-1)-roots(i)*a(2:i)
		a(1)=-roots(i)*a(1)
	end do
	d(m:1:-1)=-real(a(1:m))
	END SUBROUTINE fixrts
