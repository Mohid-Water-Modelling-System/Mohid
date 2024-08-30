	SUBROUTINE ratlsq(func,a,b,mm,kk,cof,dev)
	USE nrtype; USE nrutil, ONLY : arth,geop
	USE nr, ONLY : ratval,svbksb,svdcmp
	IMPLICIT NONE
	REAL(DP), INTENT(IN) :: a,b
	INTEGER(I4B), INTENT(IN) :: mm,kk
	REAL(DP), DIMENSION(:), INTENT(OUT) :: cof
	REAL(DP), INTENT(OUT) :: dev
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		REAL(DP), DIMENSION(:), INTENT(IN) :: x
		REAL(DP), DIMENSION(size(x)) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER(I4B), PARAMETER :: NPFAC=8,MAXIT=5
	REAL(DP), PARAMETER :: BIG=1.0e30_dp
	INTEGER(I4B) :: it,ncof,npt,npth
	REAL(DP) :: devmax,e,theta
	REAL(DP), DIMENSION((mm+kk+1)*NPFAC) :: bb,ee,fs,wt,xs
	REAL(DP), DIMENSION(mm+kk+1) :: coff,w
	REAL(DP), DIMENSION(mm+kk+1,mm+kk+1) :: v
	REAL(DP), DIMENSION((mm+kk+1)*NPFAC,mm+kk+1) :: u,temp
	ncof=mm+kk+1
	npt=NPFAC*ncof
	npth=npt/2
	dev=BIG
	theta=PIO2_D/(npt-1)
	xs(1:npth-1)=a+(b-a)*sin(theta*arth(0,1,npth-1))**2
	xs(npth:npt)=b-(b-a)*sin(theta*arth(npt-npth,-1,npt-npth+1))**2
	fs=func(xs)
	wt=1.0
	ee=1.0
	e=0.0
	do it=1,MAXIT
		bb=wt*(fs+sign(e,ee))
		temp=geop(spread(1.0_dp,1,npt),xs,ncof)
		u(:,1:mm+1)=temp(:,1:mm+1)*spread(wt,2,mm+1)
		u(:,mm+2:ncof)=-temp(:,2:ncof-mm)*spread(bb,2,ncof-mm-1)
		call svdcmp(u,w,v)
		call svbksb(u,w,v,bb,coff)
		ee=ratval(xs,coff,mm,kk)-fs
		wt=abs(ee)
		devmax=maxval(wt)
		e=sum(wt)/npt
		if (devmax <= dev) then
			cof=coff
			dev=devmax
		end if
		write(*,10) it,devmax
	end do
10	format (' ratlsq iteration=',i2,' max error=',1p,e10.3)
	END SUBROUTINE ratlsq
