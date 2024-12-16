	SUBROUTINE twofft(data1,data2,fft1,fft2)
	USE nrtype; USE nrutil, ONLY : assert,assert_eq
	USE nr, ONLY : four1
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: data1,data2
	COMPLEX(SPC), DIMENSION(:), INTENT(OUT) :: fft1,fft2
	INTEGER(I4B) :: n,n2
	COMPLEX(SPC), PARAMETER :: C1=(0.5_sp,0.0_sp), C2=(0.0_sp,-0.5_sp)
	COMPLEX, DIMENSION(size(data1)/2+1) :: h1,h2
	n=assert_eq(size(data1),size(data2),size(fft1),size(fft2),'twofft')
	call assert(iand(n,n-1)==0, 'n must be a power of 2 in twofft')
	fft1=cmplx(data1,data2,kind=spc)
	call four1(fft1,1)
	fft2(1)=cmplx(aimag(fft1(1)),0.0_sp,kind=spc)
	fft1(1)=cmplx(real(fft1(1)),0.0_sp,kind=spc)
	n2=n/2+1
	h1(2:n2)=C1*(fft1(2:n2)+conjg(fft1(n:n2:-1)))
	h2(2:n2)=C2*(fft1(2:n2)-conjg(fft1(n:n2:-1)))
	fft1(2:n2)=h1(2:n2)
	fft1(n:n2:-1)=conjg(h1(2:n2))
	fft2(2:n2)=h2(2:n2)
	fft2(n:n2:-1)=conjg(h2(2:n2))
	END SUBROUTINE twofft
