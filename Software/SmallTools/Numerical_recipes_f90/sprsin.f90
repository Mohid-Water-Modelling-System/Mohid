	SUBROUTINE sprsin_sp(a,thresh,sa)
	USE nrtype; USE nrutil, ONLY : arth,assert_eq
	IMPLICIT NONE
	REAL(SP), DIMENSION(:,:), INTENT(IN) :: a
	REAL(SP), INTENT(IN) :: thresh
	TYPE(sprs2_sp), INTENT(OUT) :: sa
	INTEGER(I4B) :: n,len
	LOGICAL(LGT), DIMENSION(size(a,1),size(a,2)) :: mask
	n=assert_eq(size(a,1),size(a,2),'sprsin_sp')
	mask=abs(a)>thresh
	len=count(mask)
	allocate(sa%val(len),sa%irow(len),sa%jcol(len))
	sa%n=n
	sa%len=len
	sa%val=pack(a,mask)
	sa%irow=pack(spread(arth(1,1,n),2,n),mask)
	sa%jcol=pack(spread(arth(1,1,n),1,n),mask)
	END SUBROUTINE sprsin_sp

	SUBROUTINE sprsin_dp(a,thresh,sa)
	USE nrtype; USE nrutil, ONLY : arth,assert_eq
	IMPLICIT NONE
	REAL(DP), DIMENSION(:,:), INTENT(IN) :: a
	REAL(DP), INTENT(IN) :: thresh
	TYPE(sprs2_dp), INTENT(OUT) :: sa
	INTEGER(I4B) :: n,len
	LOGICAL(LGT), DIMENSION(size(a,1),size(a,2)) :: mask
	n=assert_eq(size(a,1),size(a,2),'sprsin_dp')
	mask=abs(a)>thresh
	len=count(mask)
	allocate(sa%val(len),sa%irow(len),sa%jcol(len))
	sa%n=n
	sa%len=len
	sa%val=pack(a,mask)
	sa%irow=pack(spread(arth(1,1,n),2,n),mask)
	sa%jcol=pack(spread(arth(1,1,n),1,n),mask)
	END SUBROUTINE sprsin_dp
