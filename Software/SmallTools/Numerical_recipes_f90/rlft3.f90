	SUBROUTINE rlft3(data,spec,speq,isign)
	USE nrtype; USE nrutil, ONLY : assert,assert_eq
	USE nr, ONLY : four3
	REAL(SP), DIMENSION(:,:,:), INTENT(INOUT) :: data
	COMPLEX(SPC), DIMENSION(:,:,:), INTENT(INOUT) :: spec
	COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: speq
	INTEGER(I4B), INTENT(IN) :: isign
	INTEGER :: i1,i3,j1,j3,nn1,nn2,nn3
	REAL(DP) :: theta
	COMPLEX(SPC) :: c1=(0.5_sp,0.0_sp),c2,h1,h2,w
	COMPLEX(SPC), DIMENSION(size(data,2)-1) :: h1a,h2a
	COMPLEX(DPC) :: ww,wp
	c2=cmplx(0.0_sp,-0.5_sp*isign,kind=spc)
	nn1=assert_eq(size(data,1),2*size(spec,1),'rlft2: nn1')
	nn2=assert_eq(size(data,2),size(spec,2),size(speq,1),'rlft2: nn2')
	nn3=assert_eq(size(data,3),size(spec,3),size(speq,2),'rlft2: nn3')
	call assert(iand((/nn1,nn2,nn3/),(/nn1,nn2,nn3/)-1)==0, &
		'dimensions must be powers of 2 in rlft3')
	theta=TWOPI_D/(isign*nn1)
	wp=cmplx(-2.0_dp*sin(0.5_dp*theta)**2,sin(theta),kind=dpc)
	if (isign == 1) then
		spec(:,:,:)=cmplx(data(1:nn1:2,:,:),data(2:nn1:2,:,:),kind=spc)
		call four3(spec,isign)
		speq=spec(1,:,:)
	end if
	do i3=1,nn3
		j3=1
		if (i3 /= 1) j3=nn3-i3+2
		h1=c1*(spec(1,1,i3)+conjg(speq(1,j3)))
		h1a=c1*(spec(1,2:nn2,i3)+conjg(speq(nn2:2:-1,j3)))
		h2=c2*(spec(1,1,i3)-conjg(speq(1,j3)))
		h2a=c2*(spec(1,2:nn2,i3)-conjg(speq(nn2:2:-1,j3)))
		spec(1,1,i3)=h1+h2
		spec(1,2:nn2,i3)=h1a+h2a
		speq(1,j3)=conjg(h1-h2)
		speq(nn2:2:-1,j3)=conjg(h1a-h2a)
		ww=cmplx(1.0_dp,0.0_dp,kind=dpc)
		do i1=2,nn1/4+1
			j1=nn1/2-i1+2
			ww=ww*wp+ww
			w=ww
			h1=c1*(spec(i1,1,i3)+conjg(spec(j1,1,j3)))
			h1a=c1*(spec(i1,2:nn2,i3)+conjg(spec(j1,nn2:2:-1,j3)))
			h2=c2*(spec(i1,1,i3)-conjg(spec(j1,1,j3)))
			h2a=c2*(spec(i1,2:nn2,i3)-conjg(spec(j1,nn2:2:-1,j3)))
			spec(i1,1,i3)=h1+w*h2
			spec(i1,2:nn2,i3)=h1a+w*h2a
			spec(j1,1,j3)=conjg(h1-w*h2)
			spec(j1,nn2:2:-1,j3)=conjg(h1a-w*h2a)
		end do
	end do
	if (isign == -1) then
		call four3(spec,isign)
		data(1:nn1:2,:,:)=real(spec)
		data(2:nn1:2,:,:)=aimag(spec)
	end if
	END SUBROUTINE rlft3
