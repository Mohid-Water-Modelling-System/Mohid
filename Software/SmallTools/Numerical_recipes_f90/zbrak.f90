	SUBROUTINE zbrak(func,x1,x2,n,xb1,xb2,nb)
	USE nrtype; USE nrutil, ONLY : arth
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: n
	INTEGER(I4B), INTENT(OUT) :: nb
	REAL(SP), INTENT(IN) :: x1,x2
	REAL(SP), DIMENSION(:), POINTER :: xb1,xb2
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER(I4B) :: i
	REAL(SP) :: dx
	REAL(SP), DIMENSION(0:n) :: f,x
	LOGICAL(LGT), DIMENSION(1:n) :: mask
	LOGICAL(LGT), SAVE :: init=.true.
	if (init) then
		init=.false.
		nullify(xb1,xb2)
	end if
	if (associated(xb1)) deallocate(xb1)
	if (associated(xb2)) deallocate(xb2)
	dx=(x2-x1)/n
	x=x1+dx*arth(0,1,n+1)
	do i=0,n
		f(i)=func(x(i))
	end do
	mask=f(1:n)*f(0:n-1) <= 0.0
	nb=count(mask)
	allocate(xb1(nb),xb2(nb))
	xb1(1:nb)=pack(x(0:n-1),mask)
	xb2(1:nb)=pack(x(1:n),mask)
	END SUBROUTINE zbrak
