	SUBROUTINE quadvl(x,y,fa,fb,fc,fd)
	USE nrtype
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x,y
	REAL(SP), INTENT(OUT) :: fa,fb,fc,fd
	REAL(SP) :: qa,qb,qc,qd
	qa=min(2.0_sp,max(0.0_sp,1.0_sp-x))
	qb=min(2.0_sp,max(0.0_sp,1.0_sp-y))
	qc=min(2.0_sp,max(0.0_sp,x+1.0_sp))
	qd=min(2.0_sp,max(0.0_sp,y+1.0_sp))
	fa=0.25_sp*qa*qb
	fb=0.25_sp*qb*qc
	fc=0.25_sp*qc*qd
	fd=0.25_sp*qd*qa
	END SUBROUTINE quadvl
