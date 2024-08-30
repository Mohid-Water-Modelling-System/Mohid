	FUNCTION rd_s(x,y,z)
	USE nrtype; USE nrutil, ONLY : assert
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x,y,z
	REAL(SP) :: rd_s
	REAL(SP), PARAMETER :: ERRTOL=0.05_sp,TINY=1.0e-25_sp,BIG=4.5e21_sp,&
		C1=3.0_sp/14.0_sp,C2=1.0_sp/6.0_sp,C3=9.0_sp/22.0_sp,&
		C4=3.0_sp/26.0_sp,C5=0.25_sp*C3,C6=1.5_sp*C4
	REAL(SP) :: alamb,ave,delx,dely,delz,ea,eb,ec,ed,&
		ee,fac,sqrtx,sqrty,sqrtz,sum,xt,yt,zt
	call assert(min(x,y) >= 0.0, min(x+y,z) >= TINY, max(x,y,z) <= BIG, &
		'rd_s args')
	xt=x
	yt=y
	zt=z
	sum=0.0
	fac=1.0
	do
		sqrtx=sqrt(xt)
		sqrty=sqrt(yt)
		sqrtz=sqrt(zt)
		alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
		sum=sum+fac/(sqrtz*(zt+alamb))
		fac=0.25_sp*fac
		xt=0.25_sp*(xt+alamb)
		yt=0.25_sp*(yt+alamb)
		zt=0.25_sp*(zt+alamb)
		ave=0.2_sp*(xt+yt+3.0_sp*zt)
		delx=(ave-xt)/ave
		dely=(ave-yt)/ave
		delz=(ave-zt)/ave
		if (max(abs(delx),abs(dely),abs(delz)) <= ERRTOL) exit
	end do
	ea=delx*dely
	eb=delz*delz
	ec=ea-eb
	ed=ea-6.0_sp*eb
	ee=ed+ec+ec
	rd_s=3.0_sp*sum+fac*(1.0_sp+ed*(-C1+C5*ed-C6*delz*ee)&
		+delz*(C2*ee+delz*(-C3*ec+delz*C4*ea)))/(ave*sqrt(ave))
	END FUNCTION rd_s


	FUNCTION rd_v(x,y,z)
	USE nrtype; USE nrutil, ONLY : assert,assert_eq
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x,y,z
	REAL(SP), DIMENSION(size(x)) :: rd_v
	REAL(SP), PARAMETER :: ERRTOL=0.05_sp,TINY=1.0e-25_sp,BIG=4.5e21_sp,&
		C1=3.0_sp/14.0_sp,C2=1.0_sp/6.0_sp,C3=9.0_sp/22.0_sp,&
		C4=3.0_sp/26.0_sp,C5=0.25_sp*C3,C6=1.5_sp*C4
	REAL(SP), DIMENSION(size(x)) :: alamb,ave,delx,dely,delz,ea,eb,ec,ed,&
		ee,fac,sqrtx,sqrty,sqrtz,sum,xt,yt,zt
	LOGICAL(LGT), DIMENSION(size(x)) :: converged
	INTEGER(I4B) :: ndum
	ndum=assert_eq(size(x),size(y),size(z),'rd_v')
	call assert(all(min(x,y) >= 0.0), all(min(x+y,z) >= TINY), &
		all(max(x,y,z) <= BIG), 'rd_v args')
	xt=x
	yt=y
	zt=z
	sum=0.0
	fac=1.0
	converged=.false.
	do
		where (.not. converged)
			sqrtx=sqrt(xt)
			sqrty=sqrt(yt)
			sqrtz=sqrt(zt)
			alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
			sum=sum+fac/(sqrtz*(zt+alamb))
			fac=0.25_sp*fac
			xt=0.25_sp*(xt+alamb)
			yt=0.25_sp*(yt+alamb)
			zt=0.25_sp*(zt+alamb)
			ave=0.2_sp*(xt+yt+3.0_sp*zt)
			delx=(ave-xt)/ave
			dely=(ave-yt)/ave
			delz=(ave-zt)/ave
			converged = (all(max(abs(delx),abs(dely),abs(delz)) <= ERRTOL))
		end where
		if (all(converged)) exit
	end do
	ea=delx*dely
	eb=delz*delz
	ec=ea-eb
	ed=ea-6.0_sp*eb
	ee=ed+ec+ec
	rd_v=3.0_sp*sum+fac*(1.0_sp+ed*(-C1+C5*ed-C6*delz*ee)&
		+delz*(C2*ee+delz*(-C3*ec+delz*C4*ea)))/(ave*sqrt(ave))
	END FUNCTION rd_v
