	SUBROUTINE rlft2(data,spec,speq,isign)
	USE nrtype; USE nrutil, ONLY : assert,assert_eq
	USE nr, ONLY : four2
	REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: data
	COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: spec
	COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: speq
	INTEGER(I4B), INTENT(IN) :: isign
	INTEGER :: i1,j1,nn1,nn2
	REAL(DP) :: theta
	COMPLEX(SPC) :: c1=(0.5_sp,0.0_sp),c2,h1,h2,w
	COMPLEX(SPC), DIMENSION(size(data,2)-1) :: h1a,h2a
	COMPLEX(DPC) :: ww,wp
	nn1=assert_eq(size(data,1),2*size(spec,1),'rlft2: nn1')
	nn2=assert_eq(size(data,2),size(spec,2),size(speq),'rlft2: nn2')
	call assert(iand((/nn1,nn2/),(/nn1,nn2/)-1)==0, &
		'dimensions must be powers of 2 in rlft2')
	c2=cmplx(0.0_sp,-0.5_sp*isign,kind=spc)
	theta=TWOPI_D/(isign*nn1)
	wp=cmplx(-2.0_dp*sin(0.5_dp*theta)**2,sin(theta),kind=spc)
	if (isign == 1) then
		spec(:,:)=cmplx(data(1:nn1:2,:),data(2:nn1:2,:),kind=spc)
		call four2(spec,isign)
		speq=spec(1,:)
	end if
	h1=c1*(spec(1,1)+conjg(speq(1)))
	h1a=c1*(spec(1,2:nn2)+conjg(speq(nn2:2:-1)))
	h2=c2*(spec(1,1)-conjg(speq(1)))
	h2a=c2*(spec(1,2:nn2)-conjg(speq(nn2:2:-1)))
	spec(1,1)=h1+h2
	spec(1,2:nn2)=h1a+h2a
	speq(1)=conjg(h1-h2)
	speq(nn2:2:-1)=conjg(h1a-h2a)
	ww=cmplx(1.0_dp,0.0_dp,kind=dpc)
	do i1=2,nn1/4+1
		j1=nn1/2-i1+2
		ww=ww*wp+ww
		w=ww
		h1=c1*(spec(i1,1)+conjg(spec(j1,1)))
		h1a=c1*(spec(i1,2:nn2)+conjg(spec(j1,nn2:2:-1)))
		h2=c2*(spec(i1,1)-conjg(spec(j1,1)))
		h2a=c2*(spec(i1,2:nn2)-conjg(spec(j1,nn2:2:-1)))
		spec(i1,1)=h1+w*h2
		spec(i1,2:nn2)=h1a+w*h2a
		spec(j1,1)=conjg(h1-w*h2)
		spec(j1,nn2:2:-1)=conjg(h1a-w*h2a)
	end do
	if (isign == -1) then
		call four2(spec,isign)
		data(1:nn1:2,:)=real(spec)
		data(2:nn1:2,:)=aimag(spec)
	end if
	END SUBROUTINE rlft2
