	SUBROUTINE scrsho(func)
	USE nrtype
	IMPLICIT NONE
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER(I4B), PARAMETER :: ISCR=60,JSCR=21
	INTEGER(I4B) :: i,j,jz
	REAL(SP) :: dx,dyj,x,x1,x2,ybig,ysml
	REAL(SP), DIMENSION(ISCR) :: y
	CHARACTER(1), DIMENSION(ISCR,JSCR) :: scr
	CHARACTER(1) :: blank=' ',zero='-',yy='l',xx='-',ff='x'
	do
		write (*,*) ' Enter x1,x2 (= to stop)'
		read (*,*) x1,x2
		if (x1 == x2) RETURN
		scr(1,1:JSCR)=yy
		scr(ISCR,1:JSCR)=yy
		scr(2:ISCR-1,1)=xx
		scr(2:ISCR-1,JSCR)=xx
		scr(2:ISCR-1,2:JSCR-1)=blank
		dx=(x2-x1)/(ISCR-1)
		x=x1
		do i=1,ISCR
			y(i)=func(x)
			x=x+dx
		end do
		ysml=min(minval(y(:)),0.0_sp)
		ybig=max(maxval(y(:)),0.0_sp)
		if (ybig == ysml) ybig=ysml+1.0
		dyj=(JSCR-1)/(ybig-ysml)
		jz=1-ysml*dyj
		scr(1:ISCR,jz)=zero
		do i=1,ISCR
			j=1+(y(i)-ysml)*dyj
			scr(i,j)=ff
		end do
		write (*,'(1x,1p,e10.3,1x,80a1)') ybig,(scr(i,JSCR),i=1,ISCR)
		do j=JSCR-1,2,-1
			write (*,'(12x,80a1)') (scr(i,j),i=1,ISCR)
		end do
		write (*,'(1x,1p,e10.3,1x,80a1)') ysml,(scr(i,1),i=1,ISCR)
		write (*,'(12x,1p,e10.3,40x,e10.3)') x1,x2
	end do
	END SUBROUTINE scrsho
