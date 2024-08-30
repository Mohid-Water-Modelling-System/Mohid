	SUBROUTINE hufenc(ich,codep,nb,hcode)
	USE nrtype; USE nrutil, ONLY : nrerror,reallocate
	USE huf_info
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: ich
	INTEGER(I4B), INTENT(INOUT) :: nb
	CHARACTER(1), DIMENSION(:), POINTER :: codep
	TYPE(huffcode) :: hcode
	INTEGER(I4B) :: k,l,n,nc,ntmp
	k=ich+1
	if (k > hcode%nch .or. k < 1) call &
		nrerror('hufenc: ich out of range')
	do n=hcode%ncode(k),1,-1
		nc=nb/8+1
		if (nc > size(codep)) codep=>reallocate(codep,2*size(codep))
		l=mod(nb,8)
		if (l == 0) codep(nc)=char(0)
		if (btest(hcode%icode(k),n-1)) then
			ntmp=ibset(ichar(codep(nc)),l)
			codep(nc)=char(ntmp)
		end if
		nb=nb+1
	end do
	END SUBROUTINE hufenc
