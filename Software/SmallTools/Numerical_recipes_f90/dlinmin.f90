MODULE df1dim_mod
	USE nrtype
	INTEGER(I4B) :: ncom
	REAL(SP), DIMENSION(:), POINTER :: pcom,xicom
CONTAINS
!BL
	FUNCTION f1dim(x)
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: f1dim
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP) :: func
		END FUNCTION func
	END INTERFACE
	REAL(SP), DIMENSION(:), ALLOCATABLE :: xt
	allocate(xt(ncom))
	xt(:)=pcom(:)+x*xicom(:)
	f1dim=func(xt)
	deallocate(xt)
	END FUNCTION f1dim
!BL
	FUNCTION df1dim(x)
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: df1dim
	INTERFACE
		FUNCTION dfunc(x)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: dfunc
		END FUNCTION dfunc
	END INTERFACE
	REAL(SP), DIMENSION(:), ALLOCATABLE :: xt,df
	allocate(xt(ncom),df(ncom))
	xt(:)=pcom(:)+x*xicom(:)
	df(:)=dfunc(xt)
	df1dim=dot_product(df,xicom)
	deallocate(xt,df)
	END FUNCTION df1dim
END MODULE df1dim_mod

	SUBROUTINE dlinmin(p,xi,fret)
	USE nrtype; USE nrutil, ONLY : assert_eq
	USE nr, ONLY : mnbrak,dbrent
	USE df1dim_mod
	IMPLICIT NONE
	REAL(SP), INTENT(OUT) :: fret
	REAL(SP), DIMENSION(:), TARGET :: p,xi
	REAL(SP), PARAMETER :: TOL=1.0e-4_sp
	REAL(SP) :: ax,bx,fa,fb,fx,xmin,xx
	ncom=assert_eq(size(p),size(xi),'dlinmin')
	pcom=>p
	xicom=>xi
	ax=0.0
	xx=1.0
	call mnbrak(ax,xx,bx,fa,fx,fb,f1dim)
	fret=dbrent(ax,xx,bx,f1dim,df1dim,TOL,xmin)
	xi=xmin*xi
	p=p+xi
	END SUBROUTINE dlinmin
