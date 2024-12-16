	SUBROUTINE fourrow_sp(data,isign)
	USE nrtype; USE nrutil, ONLY : assert,swap
	IMPLICIT NONE
	COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: data
	INTEGER(I4B), INTENT(IN) :: isign
	INTEGER(I4B) :: n,i,istep,j,m,mmax,n2
	REAL(DP) :: theta
	COMPLEX(SPC), DIMENSION(size(data,1)) :: temp
	COMPLEX(DPC) :: w,wp
	COMPLEX(SPC) :: ws
	n=size(data,2)
	call assert(iand(n,n-1)==0, 'n must be a power of 2 in fourrow_sp')
	n2=n/2
	j=n2
	do i=1,n-2
		if (j > i) call swap(data(:,j+1),data(:,i+1))
		m=n2
		do
			if (m < 2 .or. j < m) exit
			j=j-m
			m=m/2
		end do
		j=j+m
	end do
	mmax=1
	do
		if (n <= mmax) exit
		istep=2*mmax
		theta=PI_D/(isign*mmax)
		wp=cmplx(-2.0_dp*sin(0.5_dp*theta)**2,sin(theta),kind=dpc)
		w=cmplx(1.0_dp,0.0_dp,kind=dpc)
		do m=1,mmax
			ws=w
			do i=m,n,istep
				j=i+mmax
				temp=ws*data(:,j)
				data(:,j)=data(:,i)-temp
				data(:,i)=data(:,i)+temp
			end do
			w=w*wp+w
		end do
		mmax=istep
	end do
	END SUBROUTINE fourrow_sp

	SUBROUTINE fourrow_dp(data,isign)
	USE nrtype; USE nrutil, ONLY : assert,swap
	IMPLICIT NONE
	COMPLEX(DPC), DIMENSION(:,:), INTENT(INOUT) :: data
	INTEGER(I4B), INTENT(IN) :: isign
	INTEGER(I4B) :: n,i,istep,j,m,mmax,n2
	REAL(DP) :: theta
	COMPLEX(DPC), DIMENSION(size(data,1)) :: temp
	COMPLEX(DPC) :: w,wp
	COMPLEX(DPC) :: ws
	n=size(data,2)
	call assert(iand(n,n-1)==0, 'n must be a power of 2 in fourrow_dp')
	n2=n/2
	j=n2
	do i=1,n-2
		if (j > i) call swap(data(:,j+1),data(:,i+1))
		m=n2
		do
			if (m < 2 .or. j < m) exit
			j=j-m
			m=m/2
		end do
		j=j+m
	end do
	mmax=1
	do
		if (n <= mmax) exit
		istep=2*mmax
		theta=PI_D/(isign*mmax)
		wp=cmplx(-2.0_dp*sin(0.5_dp*theta)**2,sin(theta),kind=dpc)
		w=cmplx(1.0_dp,0.0_dp,kind=dpc)
		do m=1,mmax
			ws=w
			do i=m,n,istep
				j=i+mmax
				temp=ws*data(:,j)
				data(:,j)=data(:,i)-temp
				data(:,i)=data(:,i)+temp
			end do
			w=w*wp+w
		end do
		mmax=istep
	end do
	END SUBROUTINE fourrow_dp
