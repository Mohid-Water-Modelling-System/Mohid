	SUBROUTINE zroots(a,roots,polish)
	USE nrtype; USE nrutil, ONLY : assert_eq,poly_term
	USE nr, ONLY : laguer,indexx
	IMPLICIT NONE
	COMPLEX(SPC), DIMENSION(:), INTENT(IN) :: a
	COMPLEX(SPC), DIMENSION(:), INTENT(OUT) :: roots
	LOGICAL(LGT), INTENT(IN) :: polish
	REAL(SP), PARAMETER :: EPS=1.0e-6_sp
	INTEGER(I4B) :: j,its,m
	INTEGER(I4B), DIMENSION(size(roots)) :: indx
	COMPLEX(SPC) :: x
	COMPLEX(SPC), DIMENSION(size(a)) :: ad
	m=assert_eq(size(roots),size(a)-1,'zroots')
	ad(:)=a(:)
	do j=m,1,-1
		x=cmplx(0.0_sp,kind=spc)
		call laguer(ad(1:j+1),x,its)
		if (abs(aimag(x)) <= 2.0_sp*EPS**2*abs(real(x))) &
			x=cmplx(real(x),kind=spc)
		roots(j)=x
		ad(j:1:-1)=poly_term(ad(j+1:2:-1),x)
	end do
	if (polish) then
		do j=1,m
			call laguer(a(:),roots(j),its)
		end do
	end if
	call indexx(real(roots),indx)
	roots=roots(indx)
	END SUBROUTINE zroots
