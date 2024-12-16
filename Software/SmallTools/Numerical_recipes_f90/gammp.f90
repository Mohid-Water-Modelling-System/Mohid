	FUNCTION gammp_s(a,x)
	USE nrtype; USE nrutil, ONLY : assert
	USE nr, ONLY : gcf,gser
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: a,x
	REAL(SP) :: gammp_s
	call assert( x >= 0.0,  a > 0.0, 'gammp_s args')
	if (x<a+1.0_sp) then
		gammp_s=gser(a,x)
	else
		gammp_s=1.0_sp-gcf(a,x)
	end if
	END FUNCTION gammp_s


	FUNCTION gammp_v(a,x)
	USE nrtype; USE nrutil, ONLY : assert,assert_eq
	USE nr, ONLY : gcf,gser
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: a,x
	REAL(SP), DIMENSION(size(x)) :: gammp_v
	LOGICAL(LGT), DIMENSION(size(x)) :: mask
	INTEGER(I4B) :: ndum
	ndum=assert_eq(size(a),size(x),'gammp_v')
	call assert( all(x >= 0.0),  all(a > 0.0), 'gammp_v args')
	mask = (x<a+1.0_sp)
	gammp_v=merge(gser(a,merge(x,0.0_sp,mask)), &
		1.0_sp-gcf(a,merge(x,0.0_sp,.not. mask)),mask)
	END FUNCTION gammp_v
