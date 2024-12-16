	FUNCTION ei(x)
	USE nrtype; USE nrutil, ONLY : assert,nrerror
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: ei
	INTEGER(I4B), PARAMETER :: MAXIT=100
	REAL(SP), PARAMETER :: EPS=epsilon(x),FPMIN=tiny(x)/EPS
	INTEGER(I4B) :: k
	REAL(SP) :: fact,prev,sm,term
	call assert(x > 0.0, 'ei arg')
	if (x < FPMIN) then
		ei=log(x)+EULER
	else if (x <= -log(EPS)) then
		sm=0.0
		fact=1.0
		do k=1,MAXIT
			fact=fact*x/k
			term=fact/k
			sm=sm+term
			if (term < EPS*sm) exit
		end do
		if (k > MAXIT) call nrerror('series failed in ei')
		ei=sm+log(x)+EULER
	else
		sm=0.0
		term=1.0
		do k=1,MAXIT
			prev=term
			term=term*k/x
			if (term < EPS) exit
			if (term < prev) then
				sm=sm+term
			else
				sm=sm-prev
				exit
			end if
		end do
		if (k > MAXIT) call nrerror('asymptotic failed in ei')
		ei=exp(x)*(1.0_sp+sm)/x
	end if
	END FUNCTION ei
