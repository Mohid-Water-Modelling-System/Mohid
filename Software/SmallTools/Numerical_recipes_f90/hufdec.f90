	SUBROUTINE hufdec(ich,code,nb,hcode)
	USE nrtype
	USE huf_info
	IMPLICIT NONE
	INTEGER(I4B), INTENT(OUT) :: ich
	INTEGER(I4B), INTENT(INOUT) :: nb
	CHARACTER(1), DIMENSION(:), INTENT(IN) :: code
	TYPE(huffcode) :: hcode
	INTEGER(I4B) :: l,nc,node
	node=hcode%nodemax
	do
		nc=nb/8+1
		if (nc > size(code)) then
			ich=hcode%nch
			RETURN
		end if
		l=mod(nb,8)
		nb=nb+1
		if (btest(ichar(code(nc)),l)) then
			node=hcode%iright(node)
		else
			node=hcode%left(node)
		end if
		if (node <= hcode%nch) then
			ich=node-1
			RETURN
		end if
	end do
	END SUBROUTINE hufdec
