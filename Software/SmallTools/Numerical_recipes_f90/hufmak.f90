MODULE huf_info
	USE nrtype
	IMPLICIT NONE
	TYPE huffcode
		INTEGER(I4B) :: nch,nodemax
		INTEGER(I4B), DIMENSION(:), POINTER :: icode,left,iright,ncode
	END TYPE huffcode
CONTAINS
	SUBROUTINE huff_allocate(hcode,mc)
	USE nrtype
	IMPLICIT NONE
	TYPE(huffcode) :: hcode
	INTEGER(I4B) :: mc
	INTEGER(I4B) :: mq
	mq=2*mc-1
	allocate(hcode%icode(mq),hcode%ncode(mq),hcode%left(mq),hcode%iright(mq))
	hcode%icode(:)=0
	hcode%ncode(:)=0
	END SUBROUTINE huff_allocate
!BL
	SUBROUTINE huff_deallocate(hcode)
	USE nrtype
	IMPLICIT NONE
	TYPE(huffcode) :: hcode
	deallocate(hcode%iright,hcode%left,hcode%ncode,hcode%icode)
	nullify(hcode%icode)
	nullify(hcode%ncode)
	nullify(hcode%left)
	nullify(hcode%iright)
	END SUBROUTINE huff_deallocate
END MODULE huf_info

	SUBROUTINE hufmak(nfreq,ilong,nlong,hcode)
	USE nrtype; USE nrutil, ONLY : array_copy,arth,imaxloc,nrerror
	USE huf_info
	IMPLICIT NONE
	INTEGER(I4B), INTENT(OUT) :: ilong,nlong
	INTEGER(I4B), DIMENSION(:), INTENT(IN) :: nfreq
	TYPE(huffcode) :: hcode
	INTEGER(I4B) :: ibit,j,k,n,node,nused,nerr
	INTEGER(I4B), DIMENSION(2*size(nfreq)-1) :: indx,iup,nprob
	hcode%nch=size(nfreq)
	call huff_allocate(hcode,size(nfreq))
	nused=0
	nprob(1:hcode%nch)=nfreq(1:hcode%nch)
	call array_copy(pack(arth(1,1,hcode%nch), nfreq(1:hcode%nch) /= 0 ),&
		indx,nused,nerr)
	do j=nused,1,-1
		call hufapp(j)
	end do
	k=hcode%nch
	do
		if (nused <= 1) exit
		node=indx(1)
		indx(1)=indx(nused)
		nused=nused-1
		call hufapp(1)
		k=k+1
		nprob(k)=nprob(indx(1))+nprob(node)
		hcode%left(k)=node
		hcode%iright(k)=indx(1)
		iup(indx(1))=-k
		iup(node)=k
		indx(1)=k
		call hufapp(1)
	end do
	hcode%nodemax=k
	iup(hcode%nodemax)=0
	do j=1,hcode%nch
		if (nprob(j) /= 0) then
			n=0
			ibit=0
			node=iup(j)
			do
				if (node == 0) exit
				if (node < 0) then
					n=ibset(n,ibit)
					node=-node
				end if
				node=iup(node)
				ibit=ibit+1
			end do
			hcode%icode(j)=n
			hcode%ncode(j)=ibit
		end if
	end do
	ilong=imaxloc(hcode%ncode(1:hcode%nch))
	nlong=hcode%ncode(ilong)
	if (nlong > bit_size(1_i4b)) call &
		nrerror('hufmak: Number of possible bits for code exceeded')
	CONTAINS
!BL
	SUBROUTINE hufapp(l)
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: l
	INTEGER(I4B) :: i,j,k,n
	n=nused
	i=l
	k=indx(i)
	do
		if (i > n/2) exit
		j=i+i
		if (j < n .and. nprob(indx(j)) > nprob(indx(j+1))) &
			j=j+1
		if (nprob(k) <= nprob(indx(j))) exit
		indx(i)=indx(j)
		i=j
	end do
	indx(i)=k
	END SUBROUTINE hufapp
	END SUBROUTINE hufmak
