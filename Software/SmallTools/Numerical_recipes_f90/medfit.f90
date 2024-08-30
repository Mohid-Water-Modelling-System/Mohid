	SUBROUTINE medfit(x,y,a,b,abdev)
	USE nrtype; USE nrutil, ONLY : assert_eq
	USE nr, ONLY : select
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
	REAL(SP), INTENT(OUT) :: a,b,abdev
	INTEGER(I4B) :: ndata
	REAL(SP) :: aa
	call medfit_private
	CONTAINS
!BL
	SUBROUTINE medfit_private
	IMPLICIT NONE
	REAL(SP) :: b1,b2,bb,chisq,del,f,f1,f2,sigb,sx,sxx,sxy,sy
	REAL(SP), DIMENSION(size(x)) :: tmp
	ndata=assert_eq(size(x),size(y),'medfit')
	sx=sum(x)
	sy=sum(y)
	sxy=dot_product(x,y)
	sxx=dot_product(x,x)
	del=ndata*sxx-sx**2
	aa=(sxx*sy-sx*sxy)/del
	bb=(ndata*sxy-sx*sy)/del
	tmp(:)=y(:)-(aa+bb*x(:))
	chisq=dot_product(tmp,tmp)
	sigb=sqrt(chisq/del)
	b1=bb
	f1=rofunc(b1)
	b2=bb+sign(3.0_sp*sigb,f1)
	f2=rofunc(b2)
	if (b2 == b1) then
		a=aa
		b=bb
		RETURN
	endif
	do
		if (f1*f2 <= 0.0) exit
		bb=b2+1.6_sp*(b2-b1)
		b1=b2
		f1=f2
		b2=bb
		f2=rofunc(b2)
	end do
	sigb=0.01_sp*sigb
	do
		if (abs(b2-b1) <= sigb) exit
		bb=b1+0.5_sp*(b2-b1)
		if (bb == b1 .or. bb == b2) exit
		f=rofunc(bb)
		if (f*f1 >= 0.0) then
			f1=f
			b1=bb
		else
			f2=f
			b2=bb
		end if
	end do
	a=aa
	b=bb
	abdev=abdev/ndata
	END SUBROUTINE medfit_private
!BL
	FUNCTION rofunc(b)
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: b
	REAL(SP) :: rofunc
	REAL(SP), PARAMETER :: EPS=epsilon(b)
	INTEGER(I4B) :: j
	REAL(SP), DIMENSION(size(x)) :: arr,d
	arr(:)=y(:)-b*x(:)
	if (mod(ndata,2) == 0) then
		j=ndata/2
		aa=0.5_sp*(select(j,arr)+select(j+1,arr))
	else
		aa=select((ndata+1)/2,arr)
	end if
	d(:)=y(:)-(b*x(:)+aa)
	abdev=sum(abs(d))
	where (y(:) /= 0.0) d(:)=d(:)/abs(y(:))
	rofunc=sum(x(:)*sign(1.0_sp,d(:)), mask=(abs(d(:)) > EPS) )
	END FUNCTION rofunc
	END SUBROUTINE medfit
