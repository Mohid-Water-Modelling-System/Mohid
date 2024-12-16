	FUNCTION qtrap(func,a,b)
	USE nrtype; USE nrutil, ONLY : nrerror
	USE nr, ONLY : trapzd
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: a,b
	REAL(SP) :: qtrap
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER(I4B), PARAMETER :: JMAX=20
	REAL(SP), PARAMETER :: EPS=1.0e-6_sp, UNLIKELY=-1.0e30_sp
	REAL(SP) :: olds
	INTEGER(I4B) :: j
	olds=UNLIKELY
	do j=1,JMAX
		call trapzd(func,a,b,qtrap,j)
		if (j > 5) then
			if (abs(qtrap-olds) < EPS*abs(olds) .or. &
				(qtrap == 0.0 .and. olds == 0.0)) RETURN
		end if
		olds=qtrap
	end do
	call nrerror('qtrap: too many steps')
	END FUNCTION qtrap
