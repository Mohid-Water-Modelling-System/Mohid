	SUBROUTINE four1_gather(data,isign)
	USE nrtype; USE nrutil, ONLY : arth,assert
	IMPLICIT NONE
	COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: data
	INTEGER(I4B), INTENT(IN) :: isign
	INTEGER(I4B) :: n,n2,m,mm
	REAL(DP) :: theta
	COMPLEX(SPC) :: wp
	INTEGER(I4B), DIMENSION(size(data)) :: jarr
	INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: jrev
	COMPLEX(SPC), DIMENSION(:), ALLOCATABLE :: wtab,dtemp
	n=size(data)
	call assert(iand(n,n-1)==0, 'n must be a power of 2 in four1_gather')
	if (n <= 1) RETURN
	allocate(jrev(n))
	jarr=arth(0,1,n)
	jrev=0
	n2=n/2
	m=n2
	do
		where (iand(jarr,1) /= 0) jrev=jrev+m
		jarr=jarr/2
		m=m/2
		if (m == 0) exit
	end do
	data=data(jrev+1)
	deallocate(jrev)
	allocate(dtemp(n),wtab(n2))
	jarr=arth(0,1,n)
	m=1
	mm=n2
	wtab(1)=(1.0_sp,0.0_sp)
	do
		where (iand(jarr,m) /= 0)
			dtemp=data*wtab(mm*iand(jarr,m-1)+1)
			data=eoshift(data,-m)-dtemp
		elsewhere
			data=data+eoshift(dtemp,m)
		end where
		m=m*2
		if (m >= n) exit
		mm=mm/2
		theta=PI_D/(isign*m)
		wp=cmplx(-2.0_dp*sin(0.5_dp*theta)**2, sin(theta),kind=spc)
		wtab(mm+1:n2:2*mm)=wtab(1:n2-mm:2*mm)*wp+wtab(1:n2-mm:2*mm)
	end do
	deallocate(dtemp,wtab)
	END SUBROUTINE four1_gather
