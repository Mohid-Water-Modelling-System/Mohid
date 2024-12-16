	MODULE quad3d_qgaus_mod
	USE nrtype
	PRIVATE
	PUBLIC quad3d_qgaus
	REAL(SP) :: xsav,ysav
	INTERFACE
		FUNCTION func(x,y,z)
		USE nrtype
		REAL(SP), INTENT(IN) :: x,y
		REAL(SP), DIMENSION(:), INTENT(IN) :: z
		REAL(SP), DIMENSION(size(z)) :: func
		END FUNCTION func
!BL
		FUNCTION y1(x)
		USE nrtype
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: y1
		END FUNCTION y1
!BL
		FUNCTION y2(x)
		USE nrtype
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: y2
		END FUNCTION y2
!BL
		FUNCTION z1(x,y)
		USE nrtype
		REAL(SP), INTENT(IN) :: x,y
		REAL(SP) :: z1
		END FUNCTION z1
!BL
		FUNCTION z2(x,y)
		USE nrtype
		REAL(SP), INTENT(IN) :: x,y
		REAL(SP) :: z2
		END FUNCTION z2
	END INTERFACE
	CONTAINS
!BL
	FUNCTION h(x)
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(size(x)) :: h
	INTEGER(I4B) :: i
	do i=1,size(x)
		xsav=x(i)
		h(i)=qgaus(g,y1(xsav),y2(xsav))
	end do
	END FUNCTION h
!BL
	FUNCTION g(y)
	REAL(SP), DIMENSION(:), INTENT(IN) :: y
	REAL(SP), DIMENSION(size(y)) :: g
	INTEGER(I4B) :: j
	do j=1,size(y)
		ysav=y(j)
		g(j)=qgaus(f,z1(xsav,ysav),z2(xsav,ysav))
	end do
	END FUNCTION g
!BL
	FUNCTION f(z)
	REAL(SP), DIMENSION(:), INTENT(IN) :: z
	REAL(SP), DIMENSION(size(z)) :: f
	f=func(xsav,ysav,z)
	END FUNCTION f
!BL
	RECURSIVE FUNCTION qgaus(func,a,b)
	REAL(SP), INTENT(IN) :: a,b
	REAL(SP) :: qgaus
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: func
		END FUNCTION func
	END INTERFACE
	REAL(SP) :: xm,xr
	REAL(SP), DIMENSION(5) :: dx, w = (/ 0.2955242247_sp,0.2692667193_sp,&
		0.2190863625_sp,0.1494513491_sp,0.0666713443_sp /),&
		x = (/ 0.1488743389_sp,0.4333953941_sp,0.6794095682_sp,&
		0.8650633666_sp,0.9739065285_sp /)
	xm=0.5_sp*(b+a)
	xr=0.5_sp*(b-a)
	dx(:)=xr*x(:)
	qgaus=xr*sum(w(:)*(func(xm+dx)+func(xm-dx)))
	END FUNCTION qgaus
!BL
	SUBROUTINE quad3d_qgaus(x1,x2,ss)
	REAL(SP), INTENT(IN) :: x1,x2
	REAL(SP), INTENT(OUT) :: ss
	ss=qgaus(h,x1,x2)
	END SUBROUTINE quad3d_qgaus
	END MODULE quad3d_qgaus_mod

	MODULE quad3d_qromb_mod
	USE nrtype
	PRIVATE
	PUBLIC quad3d_qromb
	REAL(SP) :: xsav,ysav
	INTERFACE
		FUNCTION func(x,y,z)
		USE nrtype
		REAL(SP), INTENT(IN) :: x,y
		REAL(SP), DIMENSION(:), INTENT(IN) :: z
		REAL(SP), DIMENSION(size(z)) :: func
		END FUNCTION func
!BL
		FUNCTION y1(x)
		USE nrtype
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: y1
		END FUNCTION y1
!BL
		FUNCTION y2(x)
		USE nrtype
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: y2
		END FUNCTION y2
!BL
		FUNCTION z1(x,y)
		USE nrtype
		REAL(SP), INTENT(IN) :: x,y
		REAL(SP) :: z1
		END FUNCTION z1
!BL
		FUNCTION z2(x,y)
		USE nrtype
		REAL(SP), INTENT(IN) :: x,y
		REAL(SP) :: z2
		END FUNCTION z2
	END INTERFACE
	CONTAINS
!BL
	FUNCTION h(x)
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(size(x)) :: h
	INTEGER(I4B) :: i
	do i=1,size(x)
		xsav=x(i)
		h(i)=qromb(g,y1(xsav),y2(xsav))
	end do
	END FUNCTION h
!BL
	FUNCTION g(y)
	REAL(SP), DIMENSION(:), INTENT(IN) :: y
	REAL(SP), DIMENSION(size(y)) :: g
	INTEGER(I4B) :: j
	do j=1,size(y)
		ysav=y(j)
		g(j)=qromb(f,z1(xsav,ysav),z2(xsav,ysav))
	end do
	END FUNCTION g
!BL
	FUNCTION f(z)
	REAL(SP), DIMENSION(:), INTENT(IN) :: z
	REAL(SP), DIMENSION(size(z)) :: f
	f=func(xsav,ysav,z)
	END FUNCTION f
!BL
	RECURSIVE FUNCTION qromb(func,a,b)
	USE nrtype; USE nrutil, ONLY : nrerror
	USE nr, ONLY : polint
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: a,b
	REAL(SP) :: qromb
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER(I4B), PARAMETER :: JMAX=20,JMAXP=JMAX+1,K=5,KM=K-1
	REAL(SP), PARAMETER :: EPS=3.0e-6_sp
	REAL(SP), DIMENSION(JMAXP) :: h,s
	REAL(SP) :: dqromb
	INTEGER(I4B) :: j
	h(1)=1.0
	do j=1,JMAX
		call trapzd(func,a,b,s(j),j)
		if (j >= K) then
			call polint(h(j-KM:j),s(j-KM:j),0.0_sp,qromb,dqromb)
			if (abs(dqromb) <= EPS*abs(qromb)) RETURN
		end if
		s(j+1)=s(j)
		h(j+1)=0.25_sp*h(j)
	end do
	call nrerror('qromb: too many steps')
	END FUNCTION qromb
!BL
	RECURSIVE SUBROUTINE trapzd(func,a,b,s,n)
	USE nrtype; USE nrutil, ONLY : arth
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: a,b
	REAL(SP), INTENT(INOUT) :: s
	INTEGER(I4B), INTENT(IN) :: n
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: func
		END FUNCTION func
	END INTERFACE
	REAL(SP) :: del,fsum
	INTEGER(I4B) :: it
	if (n == 1) then
		s=0.5_sp*(b-a)*sum(func( (/ a,b /) ))
	else
		it=2**(n-2)
		del=(b-a)/it
		fsum=sum(func(arth(a+0.5_sp*del,del,it)))
		s=0.5_sp*(s+del*fsum)
	end if
	END SUBROUTINE trapzd
!BL
	SUBROUTINE quad3d_qromb(x1,x2,ss)
	REAL(SP), INTENT(IN) :: x1,x2
	REAL(SP), INTENT(OUT) :: ss
	ss=qromb(h,x1,x2)
	END SUBROUTINE quad3d_qromb
	END MODULE quad3d_qromb_mod
