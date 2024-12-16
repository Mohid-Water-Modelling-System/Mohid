	SUBROUTINE bessjy_s(x,xnu,rj,ry,rjp,ryp)
	USE nrtype; USE nrutil, ONLY : assert,nrerror
	USE nr, ONLY : beschb
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x,xnu
	REAL(SP), INTENT(OUT) :: rj,ry,rjp,ryp
	INTEGER(I4B), PARAMETER :: MAXIT=10000
	REAL(DP), PARAMETER :: XMIN=2.0_dp,EPS=1.0e-10_dp,FPMIN=1.0e-30_dp
	INTEGER(I4B) :: i,isign,l,nl
	REAL(DP) :: a,b,c,d,del,del1,e,f,fact,fact2,fact3,ff,gam,gam1,gam2,&
		gammi,gampl,h,p,pimu,pimu2,q,r,rjl,rjl1,rjmu,rjp1,rjpl,rjtemp,&
		ry1,rymu,rymup,rytemp,sum,sum1,w,x2,xi,xi2,xmu,xmu2
	COMPLEX(DPC) :: aa,bb,cc,dd,dl,pq
	call assert(x > 0.0, xnu >= 0.0, 'bessjy args')
	nl=merge(int(xnu+0.5_dp), max(0,int(xnu-x+1.5_dp)), x < XMIN)
	xmu=xnu-nl
	xmu2=xmu*xmu
	xi=1.0_dp/x
	xi2=2.0_dp*xi
	w=xi2/PI_D
	isign=1
	h=xnu*xi
	if (h < FPMIN) h=FPMIN
	b=xi2*xnu
	d=0.0
	c=h
	do i=1,MAXIT
		b=b+xi2
		d=b-d
		if (abs(d) < FPMIN) d=FPMIN
		c=b-1.0_dp/c
		if (abs(c) < FPMIN) c=FPMIN
		d=1.0_dp/d
		del=c*d
		h=del*h
		if (d < 0.0) isign=-isign
		if (abs(del-1.0_dp) < EPS) exit
	end do
	if (i > MAXIT) call nrerror('x too large in bessjy; try asymptotic expansion')
	rjl=isign*FPMIN
	rjpl=h*rjl
	rjl1=rjl
	rjp1=rjpl
	fact=xnu*xi
	do l=nl,1,-1
		rjtemp=fact*rjl+rjpl
		fact=fact-xi
		rjpl=fact*rjtemp-rjl
		rjl=rjtemp
	end do
	if (rjl == 0.0) rjl=EPS
	f=rjpl/rjl
	if (x < XMIN) then
		x2=0.5_dp*x
		pimu=PI_D*xmu
		if (abs(pimu) < EPS) then
			fact=1.0
		else
			fact=pimu/sin(pimu)
		end if
		d=-log(x2)
		e=xmu*d
		if (abs(e) < EPS) then
			fact2=1.0
		else
			fact2=sinh(e)/e
		end if
		call beschb(xmu,gam1,gam2,gampl,gammi)
		ff=2.0_dp/PI_D*fact*(gam1*cosh(e)+gam2*fact2*d)
		e=exp(e)
		p=e/(gampl*PI_D)
		q=1.0_dp/(e*PI_D*gammi)
		pimu2=0.5_dp*pimu
		if (abs(pimu2) < EPS) then
			fact3=1.0
		else
			fact3=sin(pimu2)/pimu2
		end if
		r=PI_D*pimu2*fact3*fact3
		c=1.0
		d=-x2*x2
		sum=ff+r*q
		sum1=p
		do i=1,MAXIT
			ff=(i*ff+p+q)/(i*i-xmu2)
			c=c*d/i
			p=p/(i-xmu)
			q=q/(i+xmu)
			del=c*(ff+r*q)
			sum=sum+del
			del1=c*p-i*del
			sum1=sum1+del1
			if (abs(del) < (1.0_dp+abs(sum))*EPS) exit
		end do
		if (i > MAXIT) call nrerror('bessy series failed to converge')
		rymu=-sum
		ry1=-sum1*xi2
		rymup=xmu*xi*rymu-ry1
		rjmu=w/(rymup-f*rymu)
	else
		a=0.25_dp-xmu2
		pq=cmplx(-0.5_dp*xi,1.0_dp,kind=dpc)
		aa=cmplx(0.0_dp,xi*a,kind=dpc)
		bb=cmplx(2.0_dp*x,2.0_dp,kind=dpc)
		cc=bb+aa/pq
		dd=1.0_dp/bb
		pq=cc*dd*pq
		do i=2,MAXIT
			a=a+2*(i-1)
			bb=bb+cmplx(0.0_dp,2.0_dp,kind=dpc)
			dd=a*dd+bb
			if (absc(dd) < FPMIN) dd=FPMIN
			cc=bb+a/cc
			if (absc(cc) < FPMIN) cc=FPMIN
			dd=1.0_dp/dd
			dl=cc*dd
			pq=pq*dl
			if (absc(dl-1.0_dp) < EPS) exit
		end do
		if (i > MAXIT) call nrerror('cf2 failed in bessjy')
		p=real(pq)
		q=aimag(pq)
		gam=(p-f)/q
		rjmu=sqrt(w/((p-f)*gam+q))
		rjmu=sign(rjmu,rjl)
		rymu=rjmu*gam
		rymup=rymu*(p+q/gam)
		ry1=xmu*xi*rymu-rymup
	end if
	fact=rjmu/rjl
	rj=rjl1*fact
	rjp=rjp1*fact
	do i=1,nl
		rytemp=(xmu+i)*xi2*ry1-rymu
		rymu=ry1
		ry1=rytemp
	end do
	ry=rymu
	ryp=xnu*xi*rymu-ry1
	CONTAINS
!BL
	FUNCTION absc(z)
	IMPLICIT NONE
	COMPLEX(DPC), INTENT(IN) :: z
	REAL(DP) :: absc
	absc=abs(real(z))+abs(aimag(z))
	END FUNCTION absc
	END SUBROUTINE bessjy_s

	SUBROUTINE bessjy_v(x,xnu,rj,ry,rjp,ryp)
	USE nrtype; USE nrutil, ONLY : assert,nrerror
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: xnu
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(:), INTENT(OUT) :: rj,rjp,ry,ryp
	INTEGER(I4B), PARAMETER :: MAXIT=10000
	REAL(DP), PARAMETER :: XMIN=2.0_dp,EPS=1.0e-10_dp,FPMIN=1.0e-30_dp
	INTEGER(I4B) :: i,n
	INTEGER(I4B), DIMENSION(size(x)) :: isign
	REAL(DP), DIMENSION(size(x)) :: b,c,d,del,&
		h,xi,xi2,&
		rj_lt,ry_lt,rjp_lt,ryp_lt,rj_ge,ry_ge,rjp_ge,ryp_ge
	LOGICAL(LGT), DIMENSION(size(x)) :: converged
	call assert(all(x > 0.0), xnu >= 0.0, 'bessjy args')
	xi=1.0_dp/x
	xi2=2.0_dp*xi
	isign=1
	h=xnu*xi
	where (h < FPMIN) h=FPMIN
	b=xi2*xnu
	d=0.0
	c=h
	converged=.false.
	do i=1,MAXIT
		where (.not. converged)
			b=b+xi2
			d=b-d
			d=merge(FPMIN,d, abs(d) < FPMIN )
			c=b-1.0_dp/c
			c=merge(FPMIN,c, abs(c) < FPMIN )
			d=1.0_dp/d
			del=c*d
			h=del*h
			isign=merge(-isign,isign, d < 0.0 )
			converged=(abs(del-1.0_dp) < EPS)
		end where
		if (all(converged)) exit
	end do
	if (i > MAXIT) call nrerror('bessjy: x too large; try asymptotic expansion')
	converged=(x < XMIN)
	n=count(converged)
	call bessjy_xltxmin(pack(x,converged),xnu,pack(isign*FPMIN,converged),pack(h,converged),&
		rj_lt(1:n),ry_lt(1:n),rjp_lt(1:n),ryp_lt(1:n))
	n=size(x)-n
	call bessjy_xgexmin(pack(x,.not. converged),xnu,pack(isign,.not. converged),pack(h,.not. converged),&
		rj_ge(1:n),ry_ge(1:n),rjp_ge(1:n),ryp_ge(1:n))
	rj=unpacked(converged,rj_lt,rj_ge)
	rjp=unpacked(converged,rjp_lt,rjp_ge)
	ry=unpacked(converged,ry_lt,ry_ge)
	ryp=unpacked(converged,ryp_lt,ryp_ge)
	CONTAINS
!BL
	SUBROUTINE bessjy_xltxmin(x,xnu,rjli,hi,rj,ry,rjp,ryp)
	USE nrtype; USE nrutil, ONLY : nrerror
	USE nr, ONLY : beschb
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), INTENT(IN) :: xnu
	REAL(DP), DIMENSION(:), INTENT(IN) :: hi,rjli
	REAL(DP), DIMENSION(:), INTENT(OUT) :: rj,rjp,ry,ryp
	INTEGER(I4B) :: i,nl
	REAL(DP) :: xmu,xmu2,gam1,gam2,gammi,gampl
	REAL(DP), DIMENSION(size(x)) :: &
		c,d,del,del1,e,f,fact,fact2,fact3,ff,h,&
		p,pimu,pimu2,q,r,rjl,rjl1,rjmu,rjpl,rjp1,rjtemp,&
		ry1,rymu,rymup,rytemp,sum,sum1,w,x2,xi,xi2
	LOGICAL(LGT), DIMENSION(size(x)) :: converged
	nl=int(xnu+0.5_dp)
	xmu=xnu-nl
	xmu2=xmu*xmu
	xi=1.0_dp/x
	xi2=2.0_dp*xi
	w=xi2/PI_D
	h=hi
	where (h < FPMIN) h=FPMIN
	rjl=rjli
	rjpl=h*rjl
	rjl1=rjl
	rjp1=rjpl
	fact=xnu*xi
	do i=int(xnu+0.5_dp),1,-1
		rjtemp=fact*rjl+rjpl
		fact=fact-xi
		rjpl=fact*rjtemp-rjl
		rjl=rjtemp
	end do
	where (rjl == 0.0) rjl=EPS
	f=rjpl/rjl
	x2=0.5_dp*x
	pimu=PI_D*xmu
	where (abs(pimu) < EPS)
		fact=1.0
	elsewhere
		fact=pimu/sin(pimu)
	end where
	d=-log(x2)
	e=xmu*d
	where (abs(e) < EPS)
		fact2=1.0
	elsewhere
		fact2=sinh(e)/e
	end where
	call beschb(xmu,gam1,gam2,gampl,gammi)
	ff=2.0_dp/PI_D*fact*(gam1*cosh(e)+gam2*fact2*d)
	e=exp(e)
	p=e/(gampl*PI_D)
	q=1.0_dp/(e*PI_D*gammi)
	pimu2=0.5_dp*pimu
	where (abs(pimu2) < EPS)
		fact3=1.0
	elsewhere
		fact3=sin(pimu2)/pimu2
	end where
	r=PI_D*pimu2*fact3*fact3
	c=1.0
	d=-x2*x2
	sum=ff+r*q
	sum1=p
	converged=.false.
	do i=1,MAXIT
		where (.not. converged)
			ff=(i*ff+p+q)/(i*i-xmu2)
			c=c*d/i
			p=p/(i-xmu)
			q=q/(i+xmu)
			del=c*(ff+r*q)
			sum=sum+del
			del1=c*p-i*del
			sum1=sum1+del1
			converged=(abs(del) < (1.0_dp+abs(sum))*EPS)
		end where
		if (all(converged)) exit
	end do
	if (i > MAXIT) call nrerror('bessy series failed to converge')
	rymu=-sum
	ry1=-sum1*xi2
	rymup=xmu*xi*rymu-ry1
	rjmu=w/(rymup-f*rymu)
	do i=1,nl
		rytemp=(xmu+i)*xi2*ry1-rymu
		rymu=ry1
		ry1=rytemp
	end do
	fact=rjmu/rjl
	rj=rjl1*fact
	rjp=rjp1*fact
	ry=rymu
	ryp=xnu*xi*rymu-ry1
	END SUBROUTINE bessjy_xltxmin
!BL
	SUBROUTINE bessjy_xgexmin(x,xnu,isign,h,rj,ry,rjp,ryp)
	USE nrtype; USE nrutil, ONLY : nrerror
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), INTENT(IN) :: xnu
	INTEGER(I4B),DIMENSION(:), INTENT(IN) :: isign
	REAL(DP), DIMENSION(:), INTENT(IN) :: h
	REAL(DP), DIMENSION(:), INTENT(OUT) :: rj,rjp,ry,ryp
	INTEGER(I4B) :: i,nlmax
	INTEGER(I4B),DIMENSION(size(x)) :: nl
	REAL(DP), DIMENSION(size(x)) :: &
		a,f,fact,gam,p,q,rjl,rjl1,rjmu,rjpl,rjp1,rjtemp,ry1,rymu,rymup,rytemp,w,xi,xi2,xmu,xmu2
	COMPLEX(DPC), DIMENSION(size(x)) :: aa,bb,cc,dd,dl,pq
	LOGICAL(LGT), DIMENSION(size(x)) :: converged
	nl=max(0,int(xnu-x+1.5_dp))
	nlmax=maxval(nl)
	xmu=xnu-nl
	xi=1.0_dp/x
	xi2=2.0_dp*xi
	w=xi2/PI_D
	xmu2=xmu*xmu
	rjl=isign*FPMIN
	rjpl=h*rjl
	rjl1=rjl
	rjp1=rjpl
	fact=xnu*xi
	do i=nlmax,1,-1
		converged=(i > nl)
		if (all(converged)) exit
		where (.not. converged)
			rjtemp=fact*rjl+rjpl
			fact=fact-xi
			rjpl=fact*rjtemp-rjl
			rjl=rjtemp
		end where
	end do
	where (rjl == 0.0) rjl=EPS
	f=rjpl/rjl
	a=0.25_dp-xmu2
	pq=cmplx(-0.5_dp*xi,1.0_dp,kind=dpc)
	aa=cmplx(0.0_dp,xi*a,kind=dpc)
	bb=cmplx(2.0_dp*x,2.0_dp,kind=dpc)
	cc=bb+aa/pq
	dd=1.0_dp/bb
	pq=cc*dd*pq
	converged=.false.
	do i=2,MAXIT
		where (.not. converged)
			a=a+2*(i-1)
			bb=bb+cmplx(0.0_dp,2.0_dp,kind=dpc)
			dd=a*dd+bb
			dd=merge(cmplx(FPMIN,kind=dpc),dd, absc(dd) < FPMIN )
			cc=bb+a/cc
			cc=merge(cmplx(FPMIN,kind=dpc),cc, absc(cc) < FPMIN )
			dd=1.0_dp/dd
			dl=cc*dd
			pq=pq*dl
			converged=(absc(dl-1.0_dp) < EPS)
		end where
		if (all(converged)) exit
	end do
	if (i > MAXIT) call nrerror('bessjy: cf2 section failed')
	p=real(pq,dp)
	q=aimag(pq)
	gam=(p-f)/q
	rjmu=sqrt(w/((p-f)*gam+q))
	rjmu=sign(rjmu,rjl)
	rymu=rjmu*gam
	rymup=rymu*(p+q/gam)
	ry1=xmu*xi*rymu-rymup
	do i=1,nlmax
		converged=(i > nl)
		if (all(converged)) exit
		where (.not. converged)
			rytemp=(xmu+i)*xi2*ry1-rymu
			rymu=ry1
			ry1=rytemp
		end where
	end do
	fact=rjmu/rjl
	rj=rjl1*fact
	rjp=rjp1*fact
	ry=rymu
	ryp=xnu*xi*rymu-ry1
	END SUBROUTINE bessjy_xgexmin
!BL
	FUNCTION absc(z)
	IMPLICIT NONE
	COMPLEX(DPC), DIMENSION(:), INTENT(IN) :: z
	REAL(DP), DIMENSION(size(z)) :: absc
	absc=abs(real(z))+abs(aimag(z))
	END FUNCTION absc
!BL
	FUNCTION unpacked(mask,vtrue,vfalse)
	IMPLICIT NONE
	LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: mask
	REAL(DP), DIMENSION(:), INTENT(IN) :: vtrue, vfalse
	REAL(DP), DIMENSION(size(mask)) :: unpacked
	unpacked=unpack(vtrue,converged,0.0_dp)
	unpacked=unpack(vfalse,.not. converged,unpacked)
	END FUNCTION unpacked
	END SUBROUTINE bessjy_v
