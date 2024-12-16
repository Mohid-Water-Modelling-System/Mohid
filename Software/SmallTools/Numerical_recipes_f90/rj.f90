	FUNCTION rj_s(x,y,z,p)
	USE nrtype; USE nrutil, ONLY : assert
	USE nr, ONLY : rc,rf
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x,y,z,p
	REAL(SP) :: rj_s
	REAL(SP), PARAMETER :: ERRTOL=0.05_sp,TINY=2.5e-13_sp,BIG=9.0e11_sp,&
		C1=3.0_sp/14.0_sp,C2=1.0_sp/3.0_sp,C3=3.0_sp/22.0_sp,&
		C4=3.0_sp/26.0_sp,C5=0.75_sp*C3,C6=1.5_sp*C4,C7=0.5_sp*C2,&
		C8=C3+C3
	REAL(SP) :: a,alamb,alpha,ave,b,bet,delp,delx,&
		dely,delz,ea,eb,ec,ed,ee,fac,pt,rho,sqrtx,sqrty,sqrtz,&
		sm,tau,xt,yt,zt
	call assert(min(x,y,z) >= 0.0, min(x+y,x+z,y+z,abs(p)) >= TINY, &
		max(x,y,z,abs(p)) <= BIG, 'rj_s args')
	sm=0.0
	fac=1.0
	if (p > 0.0) then
		xt=x
		yt=y
		zt=z
		pt=p
	else
		xt=min(x,y,z)
		zt=max(x,y,z)
		yt=x+y+z-xt-zt
		a=1.0_sp/(yt-p)
		b=a*(zt-yt)*(yt-xt)
		pt=yt+b
		rho=xt*zt/yt
		tau=p*pt/yt
	end if
	do
		sqrtx=sqrt(xt)
		sqrty=sqrt(yt)
		sqrtz=sqrt(zt)
		alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
		alpha=(pt*(sqrtx+sqrty+sqrtz)+sqrtx*sqrty*sqrtz)**2
		bet=pt*(pt+alamb)**2
		sm=sm+fac*rc(alpha,bet)
		fac=0.25_sp*fac
		xt=0.25_sp*(xt+alamb)
		yt=0.25_sp*(yt+alamb)
		zt=0.25_sp*(zt+alamb)
		pt=0.25_sp*(pt+alamb)
		ave=0.2_sp*(xt+yt+zt+pt+pt)
		delx=(ave-xt)/ave
		dely=(ave-yt)/ave
		delz=(ave-zt)/ave
		delp=(ave-pt)/ave
		if (max(abs(delx),abs(dely),abs(delz),abs(delp)) <= ERRTOL) exit
	end do
	ea=delx*(dely+delz)+dely*delz
	eb=delx*dely*delz
	ec=delp**2
	ed=ea-3.0_sp*ec
	ee=eb+2.0_sp*delp*(ea-ec)
	rj_s=3.0_sp*sm+fac*(1.0_sp+ed*(-C1+C5*ed-C6*ee)+eb*(C7+delp*(-C8&
		+delp*C4))+delp*ea*(C2-delp*C3)-C2*delp*ec)/(ave*sqrt(ave))
	if (p <= 0.0) rj_s=a*(b*rj_s+3.0_sp*(rc(rho,tau)-rf(xt,yt,zt)))
	END FUNCTION rj_s


	FUNCTION rj_v(x,y,z,p)
	USE nrtype; USE nrutil, ONLY : assert,assert_eq
	USE nr, ONLY : rc,rf
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x,y,z,p
	REAL(SP), DIMENSION(size(x)) :: rj_v
	REAL(SP), PARAMETER :: ERRTOL=0.05_sp,TINY=2.5e-13_sp,BIG=9.0e11_sp,&
		C1=3.0_sp/14.0_sp,C2=1.0_sp/3.0_sp,C3=3.0_sp/22.0_sp,&
		C4=3.0_sp/26.0_sp,C5=0.75_sp*C3,C6=1.5_sp*C4,C7=0.5_sp*C2,&
		C8=C3+C3
	REAL(SP), DIMENSION(size(x)) :: a,alamb,alpha,ave,b,bet,delp,delx,&
		dely,delz,ea,eb,ec,ed,ee,fac,pt,rho,sqrtx,sqrty,sqrtz,&
		sm,tau,xt,yt,zt
	LOGICAL(LGT), DIMENSION(size(x)) :: mask
	INTEGER(I4B) :: ndum
	ndum=assert_eq(size(x),size(y),size(z),size(p),'rj_v')
	call assert(all(min(x,y,z) >= 0.0), all(min(x+y,x+z,y+z,abs(p)) >= TINY), &
		all(max(x,y,z,abs(p)) <= BIG), 'rj_v args')
	sm=0.0
	fac=1.0
	where (p > 0.0)
		xt=x
		yt=y
		zt=z
		pt=p
	elsewhere
		xt=min(x,y,z)
		zt=max(x,y,z)
		yt=x+y+z-xt-zt
		a=1.0_sp/(yt-p)
		b=a*(zt-yt)*(yt-xt)
		pt=yt+b
		rho=xt*zt/yt
		tau=p*pt/yt
	end where
	mask=.false.
	do
		where (.not. mask)
			sqrtx=sqrt(xt)
			sqrty=sqrt(yt)
			sqrtz=sqrt(zt)
			alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
			alpha=(pt*(sqrtx+sqrty+sqrtz)+sqrtx*sqrty*sqrtz)**2
			bet=pt*(pt+alamb)**2
			sm=sm+fac*rc(alpha,bet)
			fac=0.25_sp*fac
			xt=0.25_sp*(xt+alamb)
			yt=0.25_sp*(yt+alamb)
			zt=0.25_sp*(zt+alamb)
			pt=0.25_sp*(pt+alamb)
			ave=0.2_sp*(xt+yt+zt+pt+pt)
			delx=(ave-xt)/ave
			dely=(ave-yt)/ave
			delz=(ave-zt)/ave
			delp=(ave-pt)/ave
			mask = (max(abs(delx),abs(dely),abs(delz),abs(delp)) <= ERRTOL)
		end where
		if (all(mask)) exit
	end do
	ea=delx*(dely+delz)+dely*delz
	eb=delx*dely*delz
	ec=delp**2
	ed=ea-3.0_sp*ec
	ee=eb+2.0_sp*delp*(ea-ec)
	rj_v=3.0_sp*sm+fac*(1.0_sp+ed*(-C1+C5*ed-C6*ee)+eb*(C7+delp*(-C8&
		+delp*C4))+delp*ea*(C2-delp*C3)-C2*delp*ec)/(ave*sqrt(ave))
	mask = (p <= 0.0)
	where (mask) rj_v=a*(b*rj_v+&
		unpack(3.0_sp*(rc(pack(rho,mask),pack(tau,mask))-&
		rf(pack(xt,mask),pack(yt,mask),pack(zt,mask))),mask,0.0_sp))
	END FUNCTION rj_v
