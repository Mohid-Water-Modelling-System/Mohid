MODULE mpops
	USE nrtype
	INTEGER(I4B), PARAMETER :: NPAR_ICARRY=64
	CONTAINS
!BL
	SUBROUTINE icarry(karry,isum,nbits)
	IMPLICIT NONE
	INTEGER(I4B), INTENT(OUT) :: karry
	INTEGER(I2B), DIMENSION(:), INTENT(INOUT) :: isum
	INTEGER(I4B), INTENT(IN) :: nbits
	INTEGER(I4B) :: n,j
	INTEGER(I2B), DIMENSION(size(isum)) :: ihi
	INTEGER(I2B) :: mb,ihh
	n=size(isum)
	mb=ishft(1,nbits)-1
	karry=0
	if (n < NPAR_ICARRY ) then
		do j=n,2,-1
			ihh=ishft(isum(j),-nbits)
			if (ihh /= 0) then
				isum(j)=iand(isum(j),mb)
				isum(j-1)=isum(j-1)+ihh
			end if
		end do
		ihh=ishft(isum(1),-nbits)
		isum(1)=iand(isum(1),mb)
		karry=karry+ihh
	else
		do
			ihi=ishft(isum,-nbits)
			if (all(ihi == 0)) exit
			where (ihi /= 0) isum=iand(isum,mb)
			where (ihi(2:n) /= 0) isum(1:n-1)=isum(1:n-1)+ihi(2:n)
			karry=karry+ihi(1)
		end do
	end if
	END SUBROUTINE icarry
!BL
	SUBROUTINE mpadd(w,u,v,n)
	IMPLICIT NONE
	CHARACTER(1), DIMENSION(:), INTENT(OUT) :: w
	CHARACTER(1), DIMENSION(:), INTENT(IN) :: u,v
	INTEGER(I4B), INTENT(IN) :: n
	INTEGER(I2B), DIMENSION(n) :: isum
	INTEGER(I4B) :: karry
	isum=ichar(u(1:n))+ichar(v(1:n))
	call icarry(karry,isum,8_I4B)
	w(2:n+1)=char(isum)
	w(1)=char(karry)
	END SUBROUTINE mpadd
!BL
	SUBROUTINE mpsub(is,w,u,v,n)
	IMPLICIT NONE
	INTEGER(I4B), INTENT(OUT) :: is
	CHARACTER(1), DIMENSION(:), INTENT(OUT) :: w
	CHARACTER(1), DIMENSION(:), INTENT(IN) :: u,v
	INTEGER(I4B), INTENT(IN) :: n
	INTEGER(I4B) :: karry
	INTEGER(I2B), DIMENSION(n) :: isum
	isum=255+ichar(u(1:n))-ichar(v(1:n))
	isum(n)=isum(n)+1
	call icarry(karry,isum,8_I4B)
	w(1:n)=char(isum)
	is=karry-1
	END SUBROUTINE mpsub
!BL
	SUBROUTINE mpsad(w,u,n,iv)
	IMPLICIT NONE
	CHARACTER(1), DIMENSION(:), INTENT(OUT) :: w
	CHARACTER(1), DIMENSION(:), INTENT(IN) :: u
	INTEGER(I4B), INTENT(IN) :: n,iv
	INTEGER(I4B) :: karry
	INTEGER(I2B), DIMENSION(n) :: isum
	isum=ichar(u(1:n))
	isum(n)=isum(n)+iv
	call icarry(karry,isum,8_I4B)
	w(2:n+1)=char(isum)
	w(1)=char(karry)
	END SUBROUTINE mpsad
!BL
	SUBROUTINE mpsmu(w,u,n,iv)
	IMPLICIT NONE
	CHARACTER(1), DIMENSION(:), INTENT(OUT) :: w
	CHARACTER(1), DIMENSION(:), INTENT(IN) :: u
	INTEGER(I4B), INTENT(IN) :: n,iv
	INTEGER(I4B) :: karry
	INTEGER(I2B), DIMENSION(n) :: isum
	isum=ichar(u(1:n))*iv
	call icarry(karry,isum,8_I4B)
	w(2:n+1)=char(isum)
	w(1)=char(karry)
	END SUBROUTINE mpsmu
!BL
	SUBROUTINE mpneg(u,n)
	IMPLICIT NONE
	CHARACTER(1), DIMENSION(:), INTENT(INOUT) :: u
	INTEGER(I4B), INTENT(IN) :: n
	INTEGER(I4B) :: karry
	INTEGER(I2B), DIMENSION(n) :: isum
	isum=255-ichar(u(1:n))
	isum(n)=isum(n)+1
	call icarry(karry,isum,8_I4B)
	u(1:n)=char(isum)
	END SUBROUTINE mpneg
!BL
	SUBROUTINE mplsh(u,n)
	IMPLICIT NONE
	CHARACTER(1), DIMENSION(:), INTENT(INOUT) :: u
	INTEGER(I4B), INTENT(IN) :: n
	u(1:n)=u(2:n+1)
	END SUBROUTINE mplsh
!BL
	SUBROUTINE mpmov(u,v,n)
	IMPLICIT NONE
	CHARACTER(1), DIMENSION(:), INTENT(OUT) :: u
	CHARACTER(1), DIMENSION(:), INTENT(IN) :: v
	INTEGER(I4B), INTENT(IN) :: n
	u(1:n)=v(1:n)
	END SUBROUTINE mpmov
!BL
	SUBROUTINE mpsdv(w,u,n,iv,ir)
	IMPLICIT NONE
	CHARACTER(1), DIMENSION(:), INTENT(OUT) :: w
	CHARACTER(1), DIMENSION(:), INTENT(IN) :: u
	INTEGER(I4B), INTENT(IN) :: n,iv
	INTEGER(I4B), INTENT(OUT) :: ir
	INTEGER(I4B) :: i,j
	ir=0
	do j=1,n
		i=256*ir+ichar(u(j))
		w(j)=char(i/iv)
		ir=mod(i,iv)
	end do
	END SUBROUTINE mpsdv
END MODULE mpops
