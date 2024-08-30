	FUNCTION rf_s(x,y,z)
	USE nrtype; USE nrutil, ONLY : assert
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x,y,z
	REAL(SP) :: rf_s
	REAL(SP), PARAMETER :: ERRTOL=0.08_sp,TINY=1.5e-38_sp,BIG=3.0e37_sp,&
		THIRD=1.0_sp/3.0_sp,&
		C1=1.0_sp/24.0_sp,C2=0.1_sp,C3=3.0_sp/44.0_sp,C4=1.0_sp/14.0_sp
	REAL(SP) :: alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt
	call assert(min(x,y,z) >= 0.0, min(x+y,x+z,y+z) >= TINY, &
		max(x,y,z) <= BIG, 'rf_s args')
	xt=x
	yt=y
	zt=z
	do
		sqrtx=sqrt(xt)
		sqrty=sqrt(yt)
		sqrtz=sqrt(zt)
		alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
		xt=0.25_sp*(xt+alamb)
		yt=0.25_sp*(yt+alamb)
		zt=0.25_sp*(zt+alamb)
		ave=THIRD*(xt+yt+zt)
		delx=(ave-xt)/ave
		dely=(ave-yt)/ave
		delz=(ave-zt)/ave
		if (max(abs(delx),abs(dely),abs(delz)) <= ERRTOL) exit
	end do
	e2=delx*dely-delz**2
	e3=delx*dely*delz
	rf_s=(1.0_sp+(C1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(ave)
	END FUNCTION rf_s


	FUNCTION rf_v(x,y,z)
	USE nrtype; USE nrutil, ONLY : assert,assert_eq
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x,y,z
	REAL(SP), DIMENSION(size(x)) :: rf_v
	REAL(SP), PARAMETER :: ERRTOL=0.08_sp,TINY=1.5e-38_sp,BIG=3.0e37_sp,&
		THIRD=1.0_sp/3.0_sp,&
		C1=1.0_sp/24.0_sp,C2=0.1_sp,C3=3.0_sp/44.0_sp,C4=1.0_sp/14.0_sp
	REAL(SP), DIMENSION(size(x)) :: alamb,ave,delx,dely,delz,e2,e3,&
		sqrtx,sqrty,sqrtz,xt,yt,zt
	LOGICAL(LGT), DIMENSION(size(x)) :: converged
	INTEGER(I4B) :: ndum
	ndum=assert_eq(size(x),size(y),size(z),'rf_v')
	call assert(all(min(x,y,z) >= 0.0), all(min(x+y,x+z,y+z) >= TINY), &
		all(max(x,y,z) <= BIG), 'rf_v args')
	xt=x
	yt=y
	zt=z
	converged=.false.
	do
		where (.not. converged)
			sqrtx=sqrt(xt)
			sqrty=sqrt(yt)
			sqrtz=sqrt(zt)
			alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
			xt=0.25_sp*(xt+alamb)
			yt=0.25_sp*(yt+alamb)
			zt=0.25_sp*(zt+alamb)
			ave=THIRD*(xt+yt+zt)
			delx=(ave-xt)/ave
			dely=(ave-yt)/ave
			delz=(ave-zt)/ave
			converged = (max(abs(delx),abs(dely),abs(delz)) <= ERRTOL)
		end where
		if (all(converged)) exit
	end do
	e2=delx*dely-delz**2
	e3=delx*dely*delz
	rf_v=(1.0_sp+(C1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(ave)
	END FUNCTION rf_v
