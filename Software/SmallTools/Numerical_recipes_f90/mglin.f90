	SUBROUTINE mglin(u,ncycle)
	USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
	USE nr, ONLY : interp,rstrct,slvsml
	IMPLICIT NONE
	REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: u
	INTEGER(I4B), INTENT(IN) :: ncycle
	INTEGER(I4B) :: j,jcycle,n,ng,ngrid,nn
	TYPE ptr2d
		REAL(DP), POINTER :: a(:,:)
	END TYPE ptr2d
	TYPE(ptr2d), ALLOCATABLE :: rho(:)
	REAL(DP), DIMENSION(:,:), POINTER :: uj,uj_1
	n=assert_eq(size(u,1),size(u,2),'mglin')
	ng=nint(log(n-1.0)/log(2.0))
	if (n /= 2**ng+1) call nrerror('n-1 must be a power of 2 in mglin')
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
	call slvsml(uj,rho(1)%a)
	do j=2,ng
		nn=2*nn-1
		uj_1=>uj
		allocate(uj(nn,nn))
		uj=interp(uj_1)
		deallocate(uj_1)
		do jcycle=1,ncycle
			call mg(j,uj,rho(j)%a)
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
	RECURSIVE SUBROUTINE mg(j,u,rhs)
	USE nrtype
	USE nr, ONLY : interp,relax,resid,rstrct,slvsml
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: j
	REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: u
	REAL(DP), DIMENSION(:,:), INTENT(IN) :: rhs
	INTEGER(I4B), PARAMETER :: NPRE=1,NPOST=1
	INTEGER(I4B) :: jpost,jpre
	REAL(DP), DIMENSION((size(u,1)+1)/2,(size(u,1)+1)/2) :: res,v
	if (j == 1) then
		call slvsml(u,rhs)
	else
		do jpre=1,NPRE
			call relax(u,rhs)
		end do
		res=rstrct(resid(u,rhs))
		v=0.0
		call mg(j-1,v,res)
		u=u+interp(v)
		do jpost=1,NPOST
			call relax(u,rhs)
		end do
	end if
	END SUBROUTINE mg
	END SUBROUTINE mglin
