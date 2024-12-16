	SUBROUTINE jacobn(x,y,dfdx,dfdy)
	USE nrtype
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP), DIMENSION(:), INTENT(IN) :: y
	REAL(SP), DIMENSION(:), INTENT(OUT) :: dfdx
	REAL(SP), DIMENSION(:,:), INTENT(OUT) :: dfdy
	dfdx(:)=0.0
	dfdy(1,1)=-0.013_sp-1000.0_sp*y(3)
	dfdy(1,2)=0.0
	dfdy(1,3)=-1000.0_sp*y(1)
	dfdy(2,1)=0.0
	dfdy(2,2)=-2500.0_sp*y(3)
	dfdy(2,3)=-2500.0_sp*y(2)
	dfdy(3,1)=-0.013_sp-1000.0_sp*y(3)
	dfdy(3,2)=-2500.0_sp*y(3)
	dfdy(3,3)=-1000.0_sp*y(1)-2500.0_sp*y(2)
	END SUBROUTINE jacobn

	SUBROUTINE derivs(x,y,dydx)
	USE nrtype
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP), DIMENSION(:), INTENT(IN) :: y
	REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
	dydx(1)=-0.013_sp*y(1)-1000.0_sp*y(1)*y(3)
	dydx(2)=-2500.0_sp*y(2)*y(3)
	dydx(3)=-0.013_sp*y(1)-1000.0_sp*y(1)*y(3)-2500.0_sp*y(2)*y(3)
	END SUBROUTINE derivs
