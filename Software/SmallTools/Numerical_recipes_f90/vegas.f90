	SUBROUTINE vegas(region,func,init,ncall,itmx,nprn,tgral,sd,chi2a)
	USE nrtype
	USE nr, ONLY : ran1
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: region
	INTEGER(I4B), INTENT(IN) :: init,ncall,itmx,nprn
	REAL(SP), INTENT(OUT) :: tgral,sd,chi2a
	INTERFACE
		FUNCTION func(pt,wgt)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), DIMENSION(:), INTENT(IN) :: pt
		REAL(SP), INTENT(IN) :: wgt
		REAL(SP) :: func
		END FUNCTION func
	END INTERFACE
	REAL(SP), PARAMETER :: ALPH=1.5_sp,TINY=1.0e-30_sp
	INTEGER(I4B), PARAMETER :: MXDIM=10,NDMX=50
	INTEGER(I4B), SAVE :: i,it,j,k,mds,nd,ndim,ndo,ng,npg
	INTEGER(I4B), DIMENSION(MXDIM), SAVE :: ia,kg
	REAL(SP), SAVE :: calls,dv2g,dxg,f,f2,f2b,fb,rc,ti,tsi,wgt,xjac,xn,xnd,xo,harvest
	REAL(SP), DIMENSION(NDMX,MXDIM), SAVE :: d,di,xi
	REAL(SP), DIMENSION(MXDIM), SAVE :: dt,dx,x
	REAL(SP), DIMENSION(NDMX), SAVE :: r,xin
	REAL(DP), SAVE :: schi,si,swgt
	ndim=size(region)/2
	if (init <= 0) then
		mds=1
		ndo=1
		xi(1,:)=1.0
	end if
	if (init <= 1) then
		si=0.0
		swgt=0.0
		schi=0.0
	end if
	if (init <= 2) then
		nd=NDMX
		ng=1
		if (mds /= 0) then
			ng=(ncall/2.0_sp+0.25_sp)**(1.0_sp/ndim)
			mds=1
			if ((2*ng-NDMX) >= 0) then
				mds=-1
				npg=ng/NDMX+1
				nd=ng/npg
				ng=npg*nd
			end if
		end if
		k=ng**ndim
		npg=max(ncall/k,2)
		calls=real(npg,sp)*real(k,sp)
		dxg=1.0_sp/ng
		dv2g=(calls*dxg**ndim)**2/npg/npg/(npg-1.0_sp)
		xnd=nd
		dxg=dxg*xnd
		dx(1:ndim)=region(1+ndim:2*ndim)-region(1:ndim)
		xjac=1.0_sp/calls*product(dx(1:ndim))
		if (nd /= ndo) then
			r(1:max(nd,ndo))=1.0
			do j=1,ndim
				call rebin(ndo/xnd,nd,r,xin,xi(:,j))
			end do
			ndo=nd
		end if
		if (nprn >= 0) write(*,200) ndim,calls,it,itmx,nprn,&
			ALPH,mds,nd,(j,region(j),j,region(j+ndim),j=1,ndim)
	end if
	do it=1,itmx
		ti=0.0
		tsi=0.0
		kg(:)=1
		d(1:nd,:)=0.0
		di(1:nd,:)=0.0
		iterate: do
			fb=0.0
			f2b=0.0
			do k=1,npg
				wgt=xjac
				do j=1,ndim
					call ran1(harvest)
					xn=(kg(j)-harvest)*dxg+1.0_sp
					ia(j)=max(min(int(xn),NDMX),1)
					if (ia(j) > 1) then
						xo=xi(ia(j),j)-xi(ia(j)-1,j)
						rc=xi(ia(j)-1,j)+(xn-ia(j))*xo
					else
						xo=xi(ia(j),j)
						rc=(xn-ia(j))*xo
					end if
					x(j)=region(j)+rc*dx(j)
					wgt=wgt*xo*xnd
				end do
				f=wgt*func(x(1:ndim),wgt)
				f2=f*f
				fb=fb+f
				f2b=f2b+f2
				do j=1,ndim
					di(ia(j),j)=di(ia(j),j)+f
					if (mds >= 0) d(ia(j),j)=d(ia(j),j)+f2
				end do
			end do
			f2b=sqrt(f2b*npg)
			f2b=(f2b-fb)*(f2b+fb)
			if (f2b <= 0.0) f2b=TINY
			ti=ti+fb
			tsi=tsi+f2b
			if (mds < 0) then
				do j=1,ndim
					d(ia(j),j)=d(ia(j),j)+f2b
				end do
			end if
			do k=ndim,1,-1
				kg(k)=mod(kg(k),ng)+1
				if (kg(k) /= 1) cycle iterate
			end do
			exit iterate
		end do iterate
		tsi=tsi*dv2g
		wgt=1.0_sp/tsi
		si=si+real(wgt,dp)*real(ti,dp)
		schi=schi+real(wgt,dp)*real(ti,dp)**2
		swgt=swgt+real(wgt,dp)
		tgral=si/swgt
		chi2a=max((schi-si*tgral)/(it-0.99_dp),0.0_dp)
		sd=sqrt(1.0_sp/swgt)
		tsi=sqrt(tsi)
		if (nprn >= 0) then
			write(*,201) it,ti,tsi,tgral,sd,chi2a
			if (nprn /= 0) then
				do j=1,ndim
					write(*,202) j,(xi(i,j),di(i,j),&
						i=1+nprn/2,nd,nprn)
				end do
			end if
		end if
		do j=1,ndim
			xo=d(1,j)
			xn=d(2,j)
			d(1,j)=(xo+xn)/2.0_sp
			dt(j)=d(1,j)
			do i=2,nd-1
				rc=xo+xn
				xo=xn
				xn=d(i+1,j)
				d(i,j)=(rc+xn)/3.0_sp
				dt(j)=dt(j)+d(i,j)
			end do
			d(nd,j)=(xo+xn)/2.0_sp
			dt(j)=dt(j)+d(nd,j)
		end do
		where (d(1:nd,:) < TINY) d(1:nd,:)=TINY
		do j=1,ndim
			r(1:nd)=((1.0_sp-d(1:nd,j)/dt(j))/(log(dt(j))-log(d(1:nd,j))))**ALPH
			rc=sum(r(1:nd))
			call rebin(rc/xnd,nd,r,xin,xi(:,j))
		end do
	end do
200	format(/' input parameters for vegas:  ndim=',i3,'  ncall=',f8.0&
		/28x,'  it=',i5,'  itmx=',i5&
		/28x,'  nprn=',i3,'  alph=',f5.2/28x,'  mds=',i3,'   nd=',i4&
		/(30x,'xl(',i2,')= ',g11.4,' xu(',i2,')= ',g11.4))
201	format(/' iteration no.',I3,': ','integral =',g14.7,' +/- ',g9.2,&
		/' all iterations:   integral =',g14.7,' +/- ',g9.2,&
		' chi**2/it''n =',g9.2)
202	format(/' data for axis ',I2/'    X       delta i       ',&
		'   x       delta i       ','    x       delta i       ',&
		/(1x,f7.5,1x,g11.4,5x,f7.5,1x,g11.4,5x,f7.5,1x,g11.4))
	CONTAINS
!BL
	SUBROUTINE rebin(rc,nd,r,xin,xi)
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: rc
	INTEGER(I4B), INTENT(IN) :: nd
	REAL(SP), DIMENSION(:), INTENT(IN) :: r
	REAL(SP), DIMENSION(:), INTENT(OUT) :: xin
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: xi
	INTEGER(I4B) :: i,k
	REAL(SP) :: dr,xn,xo
	k=0
	xo=0.0
	dr=0.0
	do i=1,nd-1
		do
			if (rc <= dr) exit
			k=k+1
			dr=dr+r(k)
		end do
		if (k > 1) xo=xi(k-1)
		xn=xi(k)
		dr=dr-rc
		xin(i)=xn-(xn-xo)*dr/r(k)
	end do
	xi(1:nd-1)=xin(1:nd-1)
	xi(nd)=1.0
	END SUBROUTINE rebin
	END SUBROUTINE vegas
