	SUBROUTINE bandec(a,m1,m2,al,indx,d)
	USE nrtype; USE nrutil, ONLY : assert_eq,imaxloc,swap,arth
	IMPLICIT NONE
	REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
	INTEGER(I4B), INTENT(IN) :: m1,m2
	REAL(SP), DIMENSION(:,:), INTENT(OUT) :: al
	INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: indx
	REAL(SP), INTENT(OUT) :: d
	REAL(SP), PARAMETER :: TINY=1.0e-20_sp
	INTEGER(I4B) :: i,k,l,mdum,mm,n
	REAL(SP) :: dum
	n=assert_eq(size(a,1),size(al,1),size(indx),'bandec: n')
	mm=assert_eq(size(a,2),m1+m2+1,'bandec: mm')
	mdum=assert_eq(size(al,2),m1,'bandec: mdum')
	a(1:m1,:)=eoshift(a(1:m1,:),dim=2,shift=arth(m1,-1,m1))
	d=1.0
	do k=1,n
		l=min(m1+k,n)
		i=imaxloc(abs(a(k:l,1)))+k-1
		dum=a(i,1)
		if (dum == 0.0) a(k,1)=TINY
		indx(k)=i
		if (i /= k) then
			d=-d
			call swap(a(k,1:mm),a(i,1:mm))
		end if
		do i=k+1,l
			dum=a(i,1)/a(k,1)
			al(k,i-k)=dum
			a(i,1:mm-1)=a(i,2:mm)-dum*a(k,2:mm)
			a(i,mm)=0.0
		end do
	end do
	END SUBROUTINE bandec
