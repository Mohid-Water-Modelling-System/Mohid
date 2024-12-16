	SUBROUTINE bessik(x,xnu,ri,rk,rip,rkp)
	USE nrtype; USE nrutil, ONLY : assert,nrerror
	USE nr, ONLY : beschb
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x,xnu
	REAL(SP), INTENT(OUT) :: ri,rk,rip,rkp
	INTEGER(I4B), PARAMETER :: MAXIT=10000
	REAL(SP), PARAMETER :: XMIN=2.0
	REAL(DP), PARAMETER :: EPS=1.0e-10_dp,FPMIN=1.0e-30_dp
	INTEGER(I4B) :: i,l,nl
	REAL(DP) :: a,a1,b,c,d,del,del1,delh,dels,e,f,fact,fact2,ff,&
		gam1,gam2,gammi,gampl,h,p,pimu,q,q1,q2,qnew,&
		ril,ril1,rimu,rip1,ripl,ritemp,rk1,rkmu,rkmup,rktemp,&
		s,sum,sum1,x2,xi,xi2,xmu,xmu2
	call assert(x > 0.0, xnu >= 0.0, 'bessik args')
	nl=int(xnu+0.5_dp)
	xmu=xnu-nl
	xmu2=xmu*xmu
	xi=1.0_dp/x
	xi2=2.0_dp*xi
	h=xnu*xi
	if (h < FPMIN) h=FPMIN
	b=xi2*xnu
	d=0.0
	c=h
	do i=1,MAXIT
		b=b+xi2
		d=1.0_dp/(b+d)
		c=b+1.0_dp/c
		del=c*d
		h=del*h
		if (abs(del-1.0_dp) < EPS) exit
	end do
	if (i > MAXIT) call nrerror('x too large in bessik; try asymptotic expansion')
	ril=FPMIN
	ripl=h*ril
	ril1=ril
	rip1=ripl
	fact=xnu*xi
	do l=nl,1,-1
		ritemp=fact*ril+ripl
		fact=fact-xi
		ripl=fact*ritemp+ril
		ril=ritemp
	end do
	f=ripl/ril
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
		ff=fact*(gam1*cosh(e)+gam2*fact2*d)
		sum=ff
		e=exp(e)
		p=0.5_dp*e/gampl
		q=0.5_dp/(e*gammi)
		c=1.0
		d=x2*x2
		sum1=p
		do i=1,MAXIT
			ff=(i*ff+p+q)/(i*i-xmu2)
			c=c*d/i
			p=p/(i-xmu)
			q=q/(i+xmu)
			del=c*ff
			sum=sum+del
			del1=c*(p-i*ff)
			sum1=sum1+del1
			if (abs(del) < abs(sum)*EPS) exit
		end do
		if (i > MAXIT) call nrerror('bessk series failed to converge')
		rkmu=sum
		rk1=sum1*xi2
	else
		b=2.0_dp*(1.0_dp+x)
		d=1.0_dp/b
		delh=d
		h=delh
		q1=0.0
		q2=1.0
		a1=0.25_dp-xmu2
		c=a1
		q=c
		a=-a1
		s=1.0_dp+q*delh
		do i=2,MAXIT
			a=a-2*(i-1)
			c=-a*c/i
			qnew=(q1-b*q2)/a
			q1=q2
			q2=qnew
			q=q+c*qnew
			b=b+2.0_dp
			d=1.0_dp/(b+a*d)
			delh=(b*d-1.0_dp)*delh
			h=h+delh
			dels=q*delh
			s=s+dels
			if (abs(dels/s) < EPS) exit
		end do
		if (i > MAXIT) call nrerror('bessik: failure to converge in cf2')
		h=a1*h
		rkmu=sqrt(PI_D/(2.0_dp*x))*exp(-x)/s
		rk1=rkmu*(xmu+x+0.5_dp-h)*xi
	end if
	rkmup=xmu*xi*rkmu-rk1
	rimu=xi/(f*rkmu-rkmup)
	ri=(rimu*ril1)/ril
	rip=(rimu*rip1)/ril
	do i=1,nl
		rktemp=(xmu+i)*xi2*rk1+rkmu
		rkmu=rk1
		rk1=rktemp
	end do
	rk=rkmu
	rkp=xnu*xi*rkmu-rk1
	END SUBROUTINE bessik
