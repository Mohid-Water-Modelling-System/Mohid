	FUNCTION dawson_s(x)
	USE nrtype; USE nrutil, ONLY : arth,geop
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: dawson_s
	INTEGER(I4B), PARAMETER :: NMAX=6
	REAL(SP), PARAMETER :: H=0.4_sp,A1=2.0_sp/3.0_sp,A2=0.4_sp,&
		A3=2.0_sp/7.0_sp
	INTEGER(I4B) :: i,n0
	REAL(SP) :: ec,x2,xp,xx
	REAL(SP), DIMENSION(NMAX) :: d1,d2,e1
	REAL(SP), DIMENSION(NMAX), SAVE :: c=(/ (0.0_sp,i=1,NMAX) /)
	if (c(1) == 0.0) c(1:NMAX)=exp(-(arth(1,2,NMAX)*H)**2)
	if (abs(x) < 0.2_sp) then
		x2=x**2
		dawson_s=x*(1.0_sp-A1*x2*(1.0_sp-A2*x2*(1.0_sp-A3*x2)))
	else
		xx=abs(x)
		n0=2*nint(0.5_sp*xx/H)
		xp=xx-real(n0,sp)*H
		ec=exp(2.0_sp*xp*H)
		d1=arth(n0+1,2,NMAX)
		d2=arth(n0-1,-2,NMAX)
		e1=geop(ec,ec**2,NMAX)
		dawson_s=0.5641895835477563_sp*sign(exp(-xp**2),x)*&
			sum(c*(e1/d1+1.0_sp/(d2*e1)))
	end if
	END FUNCTION dawson_s

	FUNCTION dawson_v(x)
	USE nrtype; USE nrutil, ONLY : arth
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(size(x)) :: dawson_v
	INTEGER(I4B), PARAMETER :: NMAX=6
	REAL(SP), PARAMETER :: H=0.4_sp,A1=2.0_sp/3.0_sp,A2=0.4_sp,&
		A3=2.0_sp/7.0_sp
	INTEGER(I4B) :: i,n
	REAL(SP), DIMENSION(size(x)) :: x2
	REAL(SP), DIMENSION(NMAX), SAVE :: c=(/ (0.0_sp,i=1,NMAX) /)
	LOGICAL(LGT), DIMENSION(size(x)) :: mask
	if (c(1) == 0.0) c(1:NMAX)=exp(-(arth(1,2,NMAX)*H)**2)
	mask = (abs(x) >= 0.2_sp)
	dawson_v=dawsonseries_v(x,mask)
	where (.not. mask)
		x2=x**2
		dawson_v=x*(1.0_sp-A1*x2*(1.0_sp-A2*x2*(1.0_sp-A3*x2)))
	end where
	CONTAINS
!BL
	FUNCTION dawsonseries_v(xin,mask)
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: xin
	LOGICAL(LGT), DIMENSION(size(xin)), INTENT(IN) :: mask
	REAL(SP), DIMENSION(size(xin)) :: dawsonseries_v
	INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: n0
	REAL(SP), DIMENSION(:), ALLOCATABLE :: d1,d2,e1,e2,sm,xp,xx,x
	n=count(mask)
	if (n == 0) RETURN
	allocate(n0(n),d1(n),d2(n),e1(n),e2(n),sm(n),xp(n),xx(n),x(n))
	x=pack(xin,mask)
	xx=abs(x)
	n0=2*nint(0.5_sp*xx/H)
	xp=xx-real(n0,sp)*H
	e1=exp(2.0_sp*xp*H)
	e2=e1**2
	d1=n0+1.0_sp
	d2=d1-2.0_sp
	sm=0.0
	do i=1,NMAX
		sm=sm+c(i)*(e1/d1+1.0_sp/(d2*e1))
		d1=d1+2.0_sp
		d2=d2-2.0_sp
		e1=e2*e1
	end do
	sm=0.5641895835477563_sp*sign(exp(-xp**2),x)*sm
	dawsonseries_v=unpack(sm,mask,0.0_sp)
	deallocate(n0,d1,d2,e1,e2,sm,xp,xx)
	END FUNCTION dawsonseries_v
	END FUNCTION dawson_v
