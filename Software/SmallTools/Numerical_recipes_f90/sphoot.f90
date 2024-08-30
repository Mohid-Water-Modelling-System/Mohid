MODULE sphoot_data
	USE nrtype
	INTEGER(I4B) :: m,n
	REAL(SP) :: c2,dx,gamma
END MODULE sphoot_data

MODULE sphoot_caller
	USE nrtype
	INTEGER(I4B) :: nvar
	REAL(SP) :: x1,x2
END MODULE sphoot_caller

	PROGRAM sphoot
	USE nrtype; USE nrutil, ONLY : arth
	USE nr, ONLY : newt
	USE sphoot_data
	USE sphoot_caller
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NV=3,N2=1
	REAL(SP), DIMENSION(N2) :: v
	LOGICAL(LGT) :: check
	nvar=NV
	dx=1.0e-4_sp
	do
		write(*,*) 'input m,n,c-squared (999 to end)'
		read(*,*) m,n,c2
		if (c2 == 999.0) exit
		if ((n < m) .or. (m < 0)) cycle
		gamma=(-0.5_sp)**m*product(&
			arth(n+1,1,m)*(arth(real(n,sp),-1.0_sp,m)/arth(1,1,m)))
		v(1)=n*(n+1)-m*(m+1)+c2/2.0_sp
		x1=-1.0_sp+dx
		x2=0.0
		call newt(v,check)
		if (check) then
			write(*,*)'shoot failed; bad initial guess'
			exit
		else
			write(*,'(1x,t6,a)') 'mu(m,n)'
			write(*,'(1x,f12.6)') v(1)
		end if
	end do
	END PROGRAM sphoot

	SUBROUTINE load(x1,v,y)
	USE nrtype
	USE sphoot_data
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x1
	REAL(SP), DIMENSION(:), INTENT(IN) :: v
	REAL(SP), DIMENSION(:), INTENT(OUT) :: y
	REAL(SP) :: y1
	y(3)=v(1)
	y1=merge(gamma,-gamma, mod(n-m,2) == 0 )
	y(2)=-(y(3)-c2)*y1/(2*(m+1))
	y(1)=y1+y(2)*dx
	END SUBROUTINE load

	SUBROUTINE score(x2,y,f)
	USE nrtype
	USE sphoot_data
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x2
	REAL(SP), DIMENSION(:), INTENT(IN) :: y
	REAL(SP), DIMENSION(:), INTENT(OUT) :: f
	f(1)=merge(y(2),y(1), mod(n-m,2) == 0 )
	END SUBROUTINE score

	SUBROUTINE derivs(x,y,dydx)
	USE nrtype
	USE sphoot_data
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP), DIMENSION(:), INTENT(IN) :: y
	REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
	dydx(1)=y(2)
	dydx(2)=(2.0_sp*x*(m+1.0_sp)*y(2)-(y(3)-c2*x*x)*y(1))/(1.0_sp-x*x)
	dydx(3)=0.0
	END SUBROUTINE derivs
