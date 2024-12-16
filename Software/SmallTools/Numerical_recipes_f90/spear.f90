	SUBROUTINE spear(data1,data2,d,zd,probd,rs,probrs)
	USE nrtype; USE nrutil, ONLY : assert_eq
	USE nr, ONLY : betai,erfcc,sort2
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: data1,data2
	REAL(SP), INTENT(OUT) :: d,zd,probd,rs,probrs
	INTEGER(I4B) :: n
	REAL(SP) :: aved,df,en,en3n,fac,sf,sg,t,vard
	REAL(SP), DIMENSION(size(data1)) :: wksp1,wksp2
	n=assert_eq(size(data1),size(data2),'spear')
	wksp1(:)=data1(:)
	wksp2(:)=data2(:)
	call sort2(wksp1,wksp2)
	call crank(wksp1,sf)
	call sort2(wksp2,wksp1)
	call crank(wksp2,sg)
	wksp1(:)=wksp1(:)-wksp2(:)
	d=dot_product(wksp1,wksp1)
	en=n
	en3n=en**3-en
	aved=en3n/6.0_sp-(sf+sg)/12.0_sp
	fac=(1.0_sp-sf/en3n)*(1.0_sp-sg/en3n)
	vard=((en-1.0_sp)*en**2*(en+1.0_sp)**2/36.0_sp)*fac
	zd=(d-aved)/sqrt(vard)
	probd=erfcc(abs(zd)/SQRT2)
	rs=(1.0_sp-(6.0_sp/en3n)*(d+(sf+sg)/12.0_sp))/sqrt(fac)
	fac=(1.0_sp+rs)*(1.0_sp-rs)
	if (fac > 0.0) then
		t=rs*sqrt((en-2.0_sp)/fac)
		df=en-2.0_sp
		probrs=betai(0.5_sp*df,0.5_sp,df/(df+t**2))
	else
		probrs=0.0
	end if
	CONTAINS
!BL
	SUBROUTINE crank(w,s)
	USE nrtype; USE nrutil, ONLY : arth,array_copy
	IMPLICIT NONE
	REAL(SP), INTENT(OUT) :: s
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: w
	INTEGER(I4B) :: i,n,ndum,nties
	INTEGER(I4B), DIMENSION(size(w)) :: tstart,tend,tie,idx
	n=size(w)
	idx(:)=arth(1,1,n)
	tie(:)=merge(1,0,w == eoshift(w,-1))
	tie(1)=0
	w(:)=idx(:)
	if (all(tie == 0)) then
		s=0.0
		RETURN
	end if
	call array_copy(pack(idx(:),tie(:)<eoshift(tie(:),1)),tstart,nties,ndum)
	tend(1:nties)=pack(idx(:),tie(:)>eoshift(tie(:),1))
	do i=1,nties
		w(tstart(i):tend(i))=(tstart(i)+tend(i))/2.0_sp
	end do
	tend(1:nties)=tend(1:nties)-tstart(1:nties)+1
	s=sum(tend(1:nties)**3-tend(1:nties))
	END SUBROUTINE crank
	END SUBROUTINE spear
