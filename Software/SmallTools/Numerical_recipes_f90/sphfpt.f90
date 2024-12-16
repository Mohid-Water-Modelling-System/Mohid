MODULE sphfpt_data
	USE nrtype
	INTEGER(I4B) :: m,n
	REAL(SP) :: c2,dx,gamma
END MODULE sphfpt_data

MODULE sphfpt_caller
	USE nrtype
	INTEGER(I4B) :: nn2
	REAL(SP) :: x1,x2,xf
END MODULE sphfpt_caller

	PROGRAM sphfpt
	USE nrtype; USE nrutil, ONLY : arth
	USE nr, ONLY : newt
	USE sphfpt_data
	USE sphfpt_caller
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: N1=2,N2=1,NTOT=N1+N2
	REAL(SP), PARAMETER :: DXX=1.0e-4_sp
	REAL(SP), DIMENSION(:), POINTER :: v1,v2
	REAL(SP), DIMENSION(NTOT), TARGET :: v
	LOGICAL(LGT) :: check
	v1=>v(1:N2)
	v2=>v(N2+1:NTOT)
	nn2=N2
	dx=DXX
	do
		write(*,*) 'input m,n,c-squared (999 to end)'
		read(*,*) m,n,c2
		if (c2 == 999.0) exit
		if ((n < m) .or. (m < 0)) cycle
		gamma=(-0.5_sp)**m*product(&
			arth(n+1,1,m)*(arth(real(n,sp),-1.0_sp,m)/arth(1,1,m)))
		v1(1)=n*(n+1)-m*(m+1)+c2/2.0_sp
		v2(2)=v1(1)
		v2(1)=gamma*(1.0_sp-(v2(2)-c2)*dx/(2*(m+1)))
		x1=-1.0_sp+dx
		x2=1.0_sp-dx
		xf=0.0
		call newt(v,check)
		if (check) then
			write(*,*) 'shootf failed; bad initial guess'
			exit
		else
			write(*,'(1x,t6,a)') 'mu(m,n)'
			write(*,'(1x,f12.6)') v1(1)
		end if
	end do
	END PROGRAM sphfpt

	SUBROUTINE load1(x1,v1,y)
	USE nrtype
	USE sphfpt_data
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x1
	REAL(SP), DIMENSION(:), INTENT(IN) :: v1
	REAL(SP), DIMENSION(:), INTENT(OUT) :: y
	REAL(SP) :: y1
	y(3)=v1(1)
	y1=merge(gamma,-gamma,mod(n-m,2) == 0)
	y(2)=-(y(3)-c2)*y1/(2*(m+1))
	y(1)=y1+y(2)*dx
	END SUBROUTINE load1

	SUBROUTINE load2(x2,v2,y)
	USE nrtype
	USE sphfpt_data
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x2
	REAL(SP), DIMENSION(:), INTENT(IN) :: v2
	REAL(SP), DIMENSION(:), INTENT(OUT) :: y
	y(3)=v2(2)
	y(1)=v2(1)
	y(2)=(y(3)-c2)*y(1)/(2*(m+1))
	END SUBROUTINE load2

	SUBROUTINE score(xf,y,f)
	USE nrtype
	USE sphfpt_data
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: xf
	REAL(SP), DIMENSION(:), INTENT(IN) :: y
	REAL(SP), DIMENSION(:), INTENT(OUT) :: f
	f(1:3)=y(1:3)
	END SUBROUTINE score
