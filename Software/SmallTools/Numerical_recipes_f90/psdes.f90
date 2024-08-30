	SUBROUTINE psdes_s(lword,rword)
	USE nrtype
	IMPLICIT NONE
	INTEGER(I4B), INTENT(INOUT) :: lword,rword
	INTEGER(I4B), PARAMETER :: NITER=4
	INTEGER(I4B), DIMENSION(4), SAVE :: C1,C2
	DATA C1 /Z'BAA96887',Z'1E17D32C',Z'03BCDC3C',Z'0F33D1B2'/
	DATA C2 /Z'4B0F3B58',Z'E874F0C3',Z'6955C5A6',Z'55A7CA46'/
	INTEGER(I4B) :: i,ia,ib,iswap,itmph,itmpl
	do i=1,NITER
		iswap=rword
		ia=ieor(rword,C1(i))
		itmpl=iand(ia,65535)
		itmph=iand(ishft(ia,-16),65535)
		ib=itmpl**2+not(itmph**2)
		ia=ior(ishft(ib,16),iand(ishft(ib,-16),65535))
		rword=ieor(lword,ieor(C2(i),ia)+itmpl*itmph)
		lword=iswap
	end do
	END SUBROUTINE psdes_s

	SUBROUTINE psdes_v(lword,rword)
	USE nrtype; USE nrutil, ONLY : assert_eq
	IMPLICIT NONE
	INTEGER(I4B), DIMENSION(:), INTENT(INOUT) :: lword,rword
	INTEGER(I4B), PARAMETER :: NITER=4
	INTEGER(I4B), DIMENSION(4), SAVE :: C1,C2
	DATA C1 /Z'BAA96887',Z'1E17D32C',Z'03BCDC3C',Z'0F33D1B2'/
	DATA C2 /Z'4B0F3B58',Z'E874F0C3',Z'6955C5A6',Z'55A7CA46'/
	INTEGER(I4B), DIMENSION(size(lword)) :: ia,ib,iswap,itmph,itmpl
	INTEGER(I4B) :: i
	i=assert_eq(size(lword),size(rword),'psdes_v')
	do i=1,NITER
		iswap=rword
		ia=ieor(rword,C1(i))
		itmpl=iand(ia,65535)
		itmph=iand(ishft(ia,-16),65535)
		ib=itmpl**2+not(itmph**2)
		ia=ior(ishft(ib,16),iand(ishft(ib,-16),65535))
		rword=ieor(lword,ieor(C2(i),ia)+itmpl*itmph)
		lword=iswap
	end do
	END SUBROUTINE psdes_v

	SUBROUTINE psdes_safe_s(lword,rword)
	USE nrtype
	IMPLICIT NONE
	INTEGER(I4B), INTENT(INOUT) :: lword,rword
	INTEGER(I4B), PARAMETER :: NITER=4
	INTEGER(I4B), DIMENSION(4), SAVE :: C1,C2
	DATA C1 /Z'BAA96887',Z'1E17D32C',Z'03BCDC3C',Z'0F33D1B2'/
	DATA C2 /Z'4B0F3B58',Z'E874F0C3',Z'6955C5A6',Z'55A7CA46'/
	INTEGER(I4B) :: i,ia,ib,iswap
	REAL(DP) :: alo,ahi
	do i=1,NITER
		iswap=rword
		ia=ieor(rword,C1(i))
		alo=real(iand(ia,65535),dp)
		ahi=real(iand(ishft(ia,-16),65535),dp)
		ib=modint(alo*alo+real(not(modint(ahi*ahi)),dp))
		ia=ior(ishft(ib,16),iand(ishft(ib,-16),65535))
		rword=ieor(lword,modint(real(ieor(C2(i),ia),dp)+alo*ahi))
		lword=iswap
	end do
	CONTAINS
!BL
	FUNCTION modint(x)
	REAL(DP), INTENT(IN) :: x
	INTEGER(I4B) :: modint
	REAL(DP) :: a
	REAL(DP), PARAMETER :: big=huge(modint), base=big+big+2.0_dp
	a=modulo(x,base)
	if (a > big) a=a-base
	modint=nint(a,kind=i4b)
	END FUNCTION modint
	END SUBROUTINE psdes_safe_s

	SUBROUTINE psdes_safe_v(lword,rword)
	USE nrtype; USE nrutil, ONLY : assert_eq
	IMPLICIT NONE
	INTEGER(I4B), DIMENSION(:), INTENT(INOUT) :: lword,rword
	INTEGER(I4B), PARAMETER :: NITER=4
	INTEGER(I4B), SAVE :: C1(4),C2(4)
	DATA C1 /Z'BAA96887',Z'1E17D32C',Z'03BCDC3C',Z'0F33D1B2'/
	DATA C2 /Z'4B0F3B58',Z'E874F0C3',Z'6955C5A6',Z'55A7CA46'/
	INTEGER(I4B), DIMENSION(size(lword)) :: ia,ib,iswap
	REAL(DP), DIMENSION(size(lword)) :: alo,ahi
	INTEGER(I4B) :: i
	i=assert_eq(size(lword),size(rword),'psdes_safe_v')
	do i=1,NITER
		iswap=rword
		ia=ieor(rword,C1(i))
		alo=real(iand(ia,65535),dp)
		ahi=real(iand(ishft(ia,-16),65535),dp)
		ib=modint(alo*alo+real(not(modint(ahi*ahi)),dp))
		ia=ior(ishft(ib,16),iand(ishft(ib,-16),65535))
		rword=ieor(lword,modint(real(ieor(C2(i),ia),dp)+alo*ahi))
		lword=iswap
	end do
	CONTAINS
!BL
	FUNCTION modint(x)
	REAL(DP), DIMENSION(:), INTENT(IN) :: x
	INTEGER(I4B), DIMENSION(size(x)) :: modint
	REAL(DP), DIMENSION(size(x)) :: a
	REAL(DP), PARAMETER :: big=huge(modint), base=big+big+2.0_dp
	a=modulo(x,base)
	where (a > big) a=a-base
	modint=nint(a,kind=i4b)
	END FUNCTION modint
	END SUBROUTINE psdes_safe_v
