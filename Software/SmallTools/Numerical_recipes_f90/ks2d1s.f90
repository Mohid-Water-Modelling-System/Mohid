	SUBROUTINE ks2d1s(x1,y1,quadvl,d1,prob)
	USE nrtype; USE nrutil, ONLY : assert_eq
	USE nr, ONLY : pearsn,probks,quadct
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x1,y1
	REAL(SP), INTENT(OUT) :: d1,prob
	INTERFACE
		SUBROUTINE quadvl(x,y,fa,fb,fc,fd)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x,y
		REAL(SP), INTENT(OUT) :: fa,fb,fc,fd
		END SUBROUTINE quadvl
	END INTERFACE
	INTEGER(I4B) :: j,n1
	REAL(SP) :: dum,dumm,fa,fb,fc,fd,ga,gb,gc,gd,r1,rr,sqen
	n1=assert_eq(size(x1),size(y1),'ks2d1s')
	d1=0.0
	do j=1,n1
		call quadct(x1(j),y1(j),x1,y1,fa,fb,fc,fd)
		call quadvl(x1(j),y1(j),ga,gb,gc,gd)
		d1=max(d1,abs(fa-ga),abs(fb-gb),abs(fc-gc),abs(fd-gd))
	end do
	call pearsn(x1,y1,r1,dum,dumm)
	sqen=sqrt(real(n1,sp))
	rr=sqrt(1.0_sp-r1**2)
	prob=probks(d1*sqen/(1.0_sp+rr*(0.25_sp-0.75_sp/sqen)))
	END SUBROUTINE ks2d1s
