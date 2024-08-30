	SUBROUTINE spctrm(p,k,ovrlap,unit,n_window)
	USE nrtype; USE nrutil, ONLY : arth,nrerror
	USE nr, ONLY : four1
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(OUT) :: p
	INTEGER(I4B), INTENT(IN) :: k
	LOGICAL(LGT), INTENT(IN) :: ovrlap
	INTEGER(I4B), OPTIONAL, INTENT(IN) :: n_window,unit
	INTEGER(I4B) :: j,joff,joffn,kk,m,m4,m43,m44,mm,iunit,nn_window
	REAL(SP) :: den,facm,facp,sumw
	REAL(SP), DIMENSION(2*size(p)) :: w
	REAL(SP), DIMENSION(4*size(p)) :: w1
	REAL(SP), DIMENSION(size(p)) :: w2
	COMPLEX(SPC), DIMENSION(2*size(p)) :: cw1
	m=size(p)
	if (present(n_window)) then
		nn_window=n_window
	else
		nn_window=1
	end if
	if (present(unit)) then
		iunit=unit
	else
		iunit=9
	end if
	mm=m+m
	m4=mm+mm
	m44=m4+4
	m43=m4+3
	den=0.0
	facm=m
	facp=1.0_sp/m
	w1(1:mm)=window(arth(1,1,mm),facm,facp,nn_window)
	sumw=dot_product(w1(1:mm),w1(1:mm))
	p(:)=0.0
	if (ovrlap) read (iunit,*) (w2(j),j=1,m)
	do kk=1,k
		do joff=-1,0,1
			if (ovrlap) then
				w1(joff+2:joff+mm:2)=w2(1:m)
				read (iunit,*) (w2(j),j=1,m)
				joffn=joff+mm
				w1(joffn+2:joffn+mm:2)=w2(1:m)
			else
				read (iunit,*) (w1(j),j=joff+2,m4,2)
			end if
		end do
		w=window(arth(1,1,mm),facm,facp,nn_window)
		w1(2:m4:2)=w1(2:m4:2)*w
		w1(1:m4:2)=w1(1:m4:2)*w
		cw1(1:mm)=cmplx(w1(1:m4:2),w1(2:m4:2),kind=spc)
		call four1(cw1(1:mm),1)
		w1(1:m4:2)=real(cw1(1:mm))
		w1(2:m4:2)=aimag(cw1(1:mm))
		p(1)=p(1)+w1(1)**2+w1(2)**2
		p(2:m)=p(2:m)+w1(4:2*m:2)**2+w1(3:2*m-1:2)**2+&
			w1(m44-4:m44-2*m:-2)**2+w1(m43-4:m43-2*m:-2)**2
		den=den+sumw
	end do
	p(:)=p(:)/(m4*den)
	CONTAINS
!BL
	FUNCTION window(j,facm,facp,nn_window)
	IMPLICIT NONE
	INTEGER(I4B), DIMENSION(:), INTENT(IN) :: j
	INTEGER(I4B), INTENT(IN) :: nn_window
	REAL(SP), INTENT(IN) :: facm,facp
	REAL(SP), DIMENSION(size(j)) :: window
	select case(nn_window)
		case(1)
			window(j)=(1.0_sp-abs(((j-1)-facm)*facp))
		case(2)
			window(j)=1.0
		case(3)
			window(j)=(1.0_sp-(((j-1)-facm)*facp)**2)
		case default
			call nrerror('unimplemented window function in spctrm')
	end select
	END FUNCTION window
	END SUBROUTINE spctrm
