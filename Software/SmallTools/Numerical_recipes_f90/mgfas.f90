	SUBROUTINE mgfas(u,maxcyc)
	USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
	USE nr, ONLY : interp,lop,rstrct,slvsm2
	IMPLICIT NONE
	REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: u
	INTEGER(I4B), INTENT(IN) :: maxcyc
	INTEGER(I4B) :: j,jcycle,n,ng,ngrid,nn
	REAL(DP) :: res,trerr
	TYPE ptr2d
		REAL(DP), POINTER :: a(:,:)
	END TYPE ptr2d
	TYPE(ptr2d), ALLOCATABLE :: rho(:)
	REAL(DP), DIMENSION(:,:), POINTER :: uj,uj_1
	n=assert_eq(size(u,1),size(u,2),'mgfas')
	ng=nint(log(n-1.0)/log(2.0))
	if (n /= 2**ng+1) call nrerror('n-1 must be a power of 2 in mgfas')
	allocate(rho(ng))
	nn=n
	ngrid=ng
	allocate(rho(ngrid)%a(nn,nn))
	rho(ngrid)%a=u
	do
		if (nn <= 3) exit
		nn=nn/2+1
		ngrid=ngrid-1
		allocate(rho(ngrid)%a(nn,nn))
		rho(ngrid)%a=rstrct(rho(ngrid+1)%a)
	end do
	nn=3
	allocate(uj(nn,nn))
	call slvsm2(uj,rho(1)%a)
	do j=2,ng
		nn=2*nn-1
		uj_1=>uj
		allocate(uj(nn,nn))
		uj=interp(uj_1)
		deallocate(uj_1)
		do jcycle=1,maxcyc
			call mg(j,uj,trerr=trerr)
			res=sqrt(sum((lop(uj)-rho(j)%a)**2))/nn
			if (res < trerr) exit
		end do
	end do
	u=uj
	deallocate(uj)
	do j=1,ng
		deallocate(rho(j)%a)
	end do
	deallocate(rho)
	CONTAINS
!BL
	RECURSIVE SUBROUTINE mg(j,u,rhs,trerr)
	USE nrtype
	USE nr, ONLY : interp,lop,relax2,rstrct,slvsm2
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: j
	REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: u
	REAL(DP), DIMENSION(:,:), INTENT(IN), OPTIONAL :: rhs
	REAL(DP), INTENT(OUT), OPTIONAL :: trerr
	INTEGER(I4B), PARAMETER :: NPRE=1,NPOST=1
	REAL(DP), PARAMETER :: ALPHA=0.33_dp
	INTEGER(I4B) :: jpost,jpre
	REAL(DP), DIMENSION((size(u,1)+1)/2,(size(u,1)+1)/2) :: v,ut,tau
	if (j == 1) then
		call slvsm2(u,rhs+rho(j)%a)
	else
		do jpre=1,NPRE
			if (present(rhs)) then
				call relax2(u,rhs+rho(j)%a)
			else
				call relax2(u,rho(j)%a)
			end if
		end do
		ut=rstrct(u)
		v=ut
		if (present(rhs)) then
			tau=lop(ut)-rstrct(lop(u)-rhs)
		else
			tau=lop(ut)-rstrct(lop(u))
			trerr=ALPHA*sqrt(sum(tau**2))/size(tau,1)
		end if
		call mg(j-1,v,tau)
		u=u+interp(v-ut)
		do jpost=1,NPOST
			if (present(rhs)) then
				call relax2(u,rhs+rho(j)%a)
			else
				call relax2(u,rho(j)%a)
			end if
		end do
	end if
	END SUBROUTINE mg
	END SUBROUTINE mgfas
