	SUBROUTINE sort_radix(arr)
	USE nrtype; USE nrutil, ONLY : array_copy,nrerror
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
	INTEGER(I4B), DIMENSION(size(arr)) :: narr,temp
	LOGICAL, DIMENSION(size(arr)) :: msk
	INTEGER(I4B) :: k,negm,ib,ia,n,nl,nerr
	ib=bit_size(narr)
	ia=ceiling(log(real(maxexponent(arr)-minexponent(arr),sp))/log(2.0_sp)) &
		+ digits(arr)
	if (ib /= ia) call nrerror('sort_radix: bit sizes not compatible')
	negm=not(ishftc(1,-1))
	n=size(arr)
	narr=transfer(arr,narr,n)
	where (btest(narr,ib-1)) narr=ieor(narr,negm)
	do k=0,ib-2
		msk=btest(narr,k)
		call array_copy(pack(narr,.not. msk),temp,nl,nerr)
		temp(nl+1:n)=pack(narr,msk)
		narr=temp
	end do
	msk=btest(narr,ib-1)
	call array_copy(pack(narr,msk),temp,nl,nerr)
	temp(nl+1:n)=pack(narr,.not. msk)
	narr=temp
	where (btest(narr,ib-1)) narr=ieor(narr,negm)
	arr=transfer(narr,arr,n)
	END SUBROUTINE sort_radix
