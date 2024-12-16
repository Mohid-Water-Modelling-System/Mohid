	SUBROUTINE sphbes_s(n,x,sj,sy,sjp,syp)
	USE nrtype; USE nrutil, ONLY : assert
	USE nr, ONLY : bessjy
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: n
	REAL(SP), INTENT(IN) :: x
	REAL(SP), INTENT(OUT) :: sj,sy,sjp,syp
	REAL(SP), PARAMETER :: RTPIO2=1.253314137315500_sp
	REAL(SP) :: factor,order,rj,rjp,ry,ryp
	call assert(n >= 0, x > 0.0, 'sphbes_s args')
	order=n+0.5_sp
	call bessjy(x,order,rj,ry,rjp,ryp)
	factor=RTPIO2/sqrt(x)
	sj=factor*rj
	sy=factor*ry
	sjp=factor*rjp-sj/(2.0_sp*x)
	syp=factor*ryp-sy/(2.0_sp*x)
	END SUBROUTINE sphbes_s


	SUBROUTINE sphbes_v(n,x,sj,sy,sjp,syp)
	USE nrtype; USE nrutil, ONLY : assert
	USE nr, ONLY : bessjy
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: n
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(:), INTENT(OUT) :: sj,sy,sjp,syp
	REAL(SP), PARAMETER :: RTPIO2=1.253314137315500_sp
	REAL(SP) :: order
	REAL(SP), DIMENSION(size(x)) :: factor,rj,rjp,ry,ryp
	call assert(n >= 0,  all(x > 0.0), 'sphbes_v args')
	order=n+0.5_sp
	call bessjy(x,order,rj,ry,rjp,ryp)
	factor=RTPIO2/sqrt(x)
	sj=factor*rj
	sy=factor*ry
	sjp=factor*rjp-sj/(2.0_sp*x)
	syp=factor*ryp-sy/(2.0_sp*x)
	END SUBROUTINE sphbes_v
