	FUNCTION probks(alam)
	USE nrtype
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: alam
	REAL(SP) :: probks
	REAL(SP), PARAMETER :: EPS1=0.001_sp,EPS2=1.0e-8_sp
	INTEGER(I4B), PARAMETER :: NITER=100
	INTEGER(I4B) :: j
	REAL(SP) :: a2,fac,term,termbf
	a2=-2.0_sp*alam**2
	fac=2.0
	probks=0.0
	termbf=0.0
	do j=1,NITER
		term=fac*exp(a2*j**2)
		probks=probks+term
		if (abs(term) <= EPS1*termbf .or. abs(term) <= EPS2*probks) RETURN
		fac=-fac
		termbf=abs(term)
	end do
	probks=1.0
	END FUNCTION probks
