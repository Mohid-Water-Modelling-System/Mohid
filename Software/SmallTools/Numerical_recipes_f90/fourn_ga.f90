	SUBROUTINE fourn_gather(data,nn,isign)
	USE nrtype; USE nrutil, ONLY : arth,assert
	IMPLICIT NONE
	COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: data
	INTEGER(I4B), DIMENSION(:) :: nn
	INTEGER(I4B), INTENT(IN) :: isign
	INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: jarr
	INTEGER(I4B) :: ndim,idim,ntot,nprev,n,n2,msk0,msk1,msk2,m,mm,mn
	REAL(DP) :: theta
	COMPLEX(SPC) :: wp
	COMPLEX(SPC), DIMENSION(:), ALLOCATABLE :: wtab,dtemp
	call assert(iand(nn,nn-1)==0, &
		'each dimension must be a power of 2 in fourn_gather')
	ndim=size(nn)
	ntot=product(nn)
	nprev=1
	allocate(jarr(ntot))
	do idim=1,ndim
		jarr=arth(0,1,ntot)
		n=nn(idim)
		n2=n/2
		msk0=nprev
		msk1=nprev*n2
		msk2=msk0+msk1
		do
			if (msk1 <= msk0) exit
			where (iand(jarr,msk0) == 0 .neqv. iand(jarr,msk1) == 0) &
				jarr=ieor(jarr,msk2)
			msk0=msk0*2
			msk1=msk1/2
			msk2=msk0+msk1
		end do
		data=data(jarr+1)
		allocate(dtemp(ntot),wtab(n2))
		jarr=iand(n-1,arth(0,1,ntot)/nprev)
		m=1
		mm=n2
		mn=m*nprev
		wtab(1)=(1.0_sp,0.0_sp)
		do
			if (mm == 0) exit
			where (iand(jarr,m) /= 0)
				dtemp=data*wtab(mm*iand(jarr,m-1)+1)
				data=eoshift(data,-mn)-dtemp
			elsewhere
				data=data+eoshift(dtemp,mn)
			end where
			m=m*2
			if (m >= n) exit
			mn=m*nprev
			mm=mm/2
			theta=PI_D/(isign*m)
			wp=cmplx(-2.0_dp*sin(0.5_dp*theta)**2,sin(theta),kind=dpc)
			wtab(mm+1:n2:2*mm)=wtab(1:n2-mm:2*mm)*wp &
				+wtab(1:n2-mm:2*mm)
		end do
		deallocate(dtemp,wtab)
		nprev=n*nprev
	end do
	deallocate(jarr)
	END SUBROUTINE fourn_gather
