	FUNCTION icrc(crc,buf,jinit,jrev)
	USE nrtype
	IMPLICIT NONE
	CHARACTER(1), DIMENSION(:), INTENT(IN) :: buf
	INTEGER(I2B), INTENT(IN) :: crc,jinit
	INTEGER(I4B), INTENT(IN) :: jrev
	INTEGER(I2B) :: icrc
	INTEGER(I4B), SAVE :: init=0
	INTEGER(I2B) :: j,cword,ich
	INTEGER(I2B), DIMENSION(0:255), SAVE :: icrctb,rchr
	INTEGER(I2B), DIMENSION(0:15) :: it = &
		(/ 0,8,4,12,2,10,6,14,1,9,5,13,3,11,7,15 /)
	if (init == 0) then
		init=1
		do j=0,255
			icrctb(j)=icrc1(ishft(j,8),char(0))
			rchr(j)=ishft(it(iand(j,15_I2B)),4)+it(ishft(j,-4))
		end do
	end if
	cword=crc
	if (jinit >= 0) then
		cword=ior(jinit,ishft(jinit,8))
	else if (jrev < 0) then
		cword=ior(rchr(hibyte()),ishft(rchr(lobyte()),8))
	end if
	do j=1,size(buf)
		ich=ichar(buf(j))
		if (jrev < 0) ich=rchr(ich)
		cword=ieor(icrctb(ieor(ich,hibyte())),ishft(lobyte(),8))
	end do
	icrc=merge(cword, &
		ior(rchr(hibyte()),ishft(rchr(lobyte()),8)), jrev >= 0)
	CONTAINS
!BL
	FUNCTION hibyte()
	INTEGER(I2B) :: hibyte
	hibyte = ishft(cword,-8)
	END FUNCTION hibyte
!BL
	FUNCTION lobyte()
	INTEGER(I2B) :: lobyte
	lobyte = iand(cword,255_I2B)
	END FUNCTION lobyte
!BL
	FUNCTION icrc1(crc,onech)
	INTEGER(I2B), INTENT(IN) :: crc
	CHARACTER(1), INTENT(IN) :: onech
	INTEGER(I2B) :: icrc1
	INTEGER(I2B) :: i,ich, bit16, ccitt
	DATA bit16,ccitt /Z'8000', Z'1021'/
	ich=ichar(onech)
	icrc1=ieor(crc,ishft(ich,8))
	do i=1,8
		icrc1=merge(ieor(ccitt,ishft(icrc1,1)), &
			ishft(icrc1,1), iand(icrc1,bit16) /= 0)
	end do
	END FUNCTION icrc1
	END FUNCTION icrc
